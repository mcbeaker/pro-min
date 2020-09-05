
class SlowBar(Bar):
    suffix = '%(index)d/%(max)d - %(remaining_hours)d hours remaining'
    @property
    def remaining_hours(self):
        return self.eta // 3600