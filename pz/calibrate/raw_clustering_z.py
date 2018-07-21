
from astropy.cosmology import Plank15
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as inter_spline
from scipy.spatial import cKDTree

from ceci import PipelineStage
from descformats import \
    HDFFile, MetacalCatalog, TomographyCatalog, RandomsCatalog, YamlFile


"""Compute the number over-density of galaxies from a tomographic bin around
"""


class RawClusteringZ(PipelineStage):
    name = "RawClusteringZ"
    # Need to load multiple
    inputs = [
        ('random_catalog', RandomsCatalog),
        ('reference_catalog', HDFFile)
        ('shear_catalog', MetacalCatalog),
        ('tomogrphy_catalog', TomographyCatalog),
    ]
    outputs = [
        ('raw_clustering_z', HDFFile),
    ]
    required_config = {
        'min_sep_mpc': 0.1,
        'max_sep_mpc': 1.0,
    }

    def __init__(**kwargs):
        PipelineStage.__init__(**kwargs)
        self._create_angle_dist_splines()

    def run(self):

        ref_points, ref_zs = self.load_reference_catalog()
        unkn_points, unkn_weights = self.load_unknown_catalog()
        rand_points = self.load_random_catalog()
        rand_weights = unkn_weights[np.random.randint(len(unkn_weights),
                                                      size=len(rand_points))]

        unkn_tree = cKDTree(unkn_points)
        rand_tree = cKDTree(ref_points)

        ref_n_unkn_pair = np.zeros_like(ref_zs)
        ref_n_rand_pair = np.zeros_like(ref_zs)

        for ref_idx, (ref_point, ref_z) in enumerate(zip(ref_points, ref_zs)):

            comov_dist = Plank15.comoving_distance(ref_z).value
            min_sep_cos = np.cos(self.min_sep_mpc / comov_dist)
            max_sep_dist = self._angle_to_dist(self.max_sep_mpc / comov_dist)

            ref_n_unkn_pair[ref_idx] = self._compute_dist_weighted_pairs(
                ref_point, unkn_tree, unkn_points, unkn_weights, comov_dist,
                min_sep_cos, max_sep_dist)
            ref_n_rand_pair[ref_idx] = self._compute_dist_weighted_pairs(
                ref_point, unkn_tree, unkn_points, unkn_weights, comov_dist,
                min_sep_cos, max_sep_dist)


    def _compute_dist_weighted_pairs(self, point, tree, tree_points, weights,
                                     comov_dist, min_sep_cos, max_sep_dist):

        near_pair_ids = np.array(tree.query_ball_point(point, max_sep_dist))
        pair_cosines = np.dot(tree_points[near_pair_ids], point)

        args = np.argwhere(pair_cosines > min_sep_cos)
        near_pair_ids = near_pair_ids[args]
        pair_mpcs = np.arccos(pair_cosines[args]) * comov_dist

        return np.sum(weights[near_pair_ids] / pair_mpcs)

    def _create_angle_spline(self):
        angles = np.linspace(0, np.pi / 2, 10000)
        dists = np.sqrt(2 - 2 * np.cos(angles))
        self._dist_to_angle = inter_spline(angles, dists)
        self._angle_to_dist = inter_spline(dists, angles)

    def load_reference_catalog(self):
        pass

    def load_unknown_catalog(self):
        pass

    def load_random_catalog(self):
        pass


if __name__ == '__main__':
    PipelineStage.main()
