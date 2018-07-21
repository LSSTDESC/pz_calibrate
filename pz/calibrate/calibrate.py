
from ceci import PipelineStage
from descformats import \
    TextFile, FitsFile, HDFFile, YamlFile, TomographyCatalog, RandomsCatalog, \
    HDFFile


"""Combine raw N(z)
"""


class PZCalibrate(PipelineStage):
    name = "PZCalibrate"
    # Need to load multiple
    inputs = [
        ('tomographic_bin1', YamlFile),
        ('some_input_tag', TextFile),
    ]
    outputs = [
        ('number_vs_redshift', HDFFile),
        # More inputs can go here
    ]
    required_config = {
    }

    def run(self):
        pass
