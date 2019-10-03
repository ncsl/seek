class Volume(object):
    origin = []
    voxel_size = -1
    voxel_unit = ""

    def _find_summary_info(self):
        summary = {
            "Volume type": self.__class__.__name__,
            "Origin": self.origin,
            "Voxel size": self.voxel_size,
            "Units": self.voxel_unit,
        }
        return summary
