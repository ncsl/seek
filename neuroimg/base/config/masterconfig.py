from neuroimg.base.config.config import GenericConfig, FiguresConfig, \
    CalculusConfig, InputConfig, OutputConfig


class Config(object):
    generic = GenericConfig()
    figures = FiguresConfig()
    calcul = CalculusConfig()

    def __init__(self,
                 raw_data_folder=None,
                 output_base=None,
                 separate_by_run=False):
        self.input = InputConfig(raw_data_folder)
        self.out = OutputConfig(output_base, separate_by_run)
