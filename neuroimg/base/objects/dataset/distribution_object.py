from edp.base.objects.dataset.basedataobject import BaseDataset


class Distribution(BaseDataset):
    """
    A class for storing distribution results from multiple datasets.

    This can be for multiple patients, or multiple datasets.
    """

    def __init__(self, distributions, label_groups, outcome_list=[], engelscore_list=[]):
        self.distributions = distributions
        self.label_groups = label_groups
        self.outcome_list = outcome_list
        self.engelscore_list = engelscore_list

    @property
    def groups(self):
        return self.label_groups.keys()

    @property
    def label_inds(self):
        return self.label_groups.values()

    @property
    def success_inds(self):
        return [idx for idx, outcome in enumerate(self.outcome_list) if outcome == 's']

    @property
    def fail_inds(self):
        return [idx for idx, outcome in enumerate(self.outcome_list) if outcome == 'f']

    @property
    def engel1_inds(self):
        return [idx for idx, score in enumerate(self.engelscore_list) if score == 1]

    def get_engelinds(self, engelscore):
        return [idx for idx, score in enumerate(self.engelscore_list) if score == engelscore]
