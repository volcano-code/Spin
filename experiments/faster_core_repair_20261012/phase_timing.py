"""Diagnostic-only timings for the existing native stage observer."""
import time


class PhaseTiming:
    def __init__(self, clock=time.perf_counter):
        self.clock = clock
        self.active = {}
        self.stages = {}
        self.events = []

    def __call__(self, name):
        now = self.clock()
        self.events.append(dict(stage=name, seconds=now))
        if name.endswith('-begin'):
            key = name[:-6]
            if key in self.active:
                raise ValueError('DUPLICATE_STAGE_BEGIN:'+key)
            self.active[key] = now
        elif name.endswith('-end'):
            key = name[:-4]
            if key not in self.active:
                raise ValueError('STAGE_END_WITHOUT_BEGIN:'+key)
            elapsed = now-self.active.pop(key)
            if elapsed < 0:
                raise ValueError('NON_MONOTONIC_STAGE_CLOCK')
            record = self.stages.setdefault(key, dict(calls=0, seconds=0.))
            record['calls'] += 1
            record['seconds'] += elapsed

    def finish(self):
        return dict(stages=self.stages, incomplete_stages=sorted(self.active), events=self.events,
                    interpretation='nested stage times must not be summed as disjoint costs')
