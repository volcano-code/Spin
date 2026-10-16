"""Exclusive formal records with explicit durability and measured I/O costs."""
import json
import os
from pathlib import Path
import time


class RunJournal:
    def __init__(self,root,*,durability_profile='PER_EVENT_FSYNC_V1',json_value=None):
        if durability_profile not in ('PER_EVENT_FSYNC_V1','EVALUATION_CHECKPOINT_V1'):
            raise ValueError('unknown run durability profile')
        self.root=Path(root)
        if not self.root.is_dir() or self.root.is_symlink():
            raise ValueError('run journal requires an existing real directory')
        self.durability_profile=durability_profile
        self.json_value=json_value
        self.metrics={'serialization_ns':0,'write_flush_ns':0,
            'event_fsync_ns':0,'event_fsync_calls':0,'file_fsync_ns':0,'file_fsync_calls':0,
            'directory_fsync_ns':0,'directory_fsync_calls':0}
        self._pending_directory=False
        self._closed=False
        self.events=(self.root/'events.jsonl').open('x')

    def __enter__(self):return self

    def __exit__(self,kind,value,traceback):
        self.close()
        return False

    def _require_open(self):
        if self._closed:raise ValueError('run journal is closed')

    def _fsync(self,fd,kind):
        start=time.monotonic_ns()
        try:os.fsync(fd)
        finally:
            self.metrics[kind+'_fsync_ns']+=time.monotonic_ns()-start
            self.metrics[kind+'_fsync_calls']+=1

    def _directory_sync(self):
        descriptor=os.open(self.root,os.O_RDONLY|os.O_DIRECTORY)
        try:self._fsync(descriptor,'directory')
        finally:os.close(descriptor)
        self._pending_directory=False

    def _serialize(self,value):
        start=time.monotonic_ns()
        try:return json.dumps(value,default=self.json_value,sort_keys=True,allow_nan=False)
        finally:self.metrics['serialization_ns']+=time.monotonic_ns()-start

    def publish(self,name,value):
        self._require_open()
        if not isinstance(name,str) or not name or Path(name).name!=name:
            raise ValueError('run record name must be a single filename')
        payload=self._serialize(value)
        with (self.root/name).open('x') as stream:
            start=time.monotonic_ns();stream.write(payload);stream.flush()
            self.metrics['write_flush_ns']+=time.monotonic_ns()-start
            self._fsync(stream.fileno(),'file')
        if self.durability_profile=='EVALUATION_CHECKPOINT_V1' and name.startswith('evaluation-'):
            self._pending_directory=True
        else:self._directory_sync()

    def emit(self,row):
        self._require_open()
        if not isinstance(row,dict):raise ValueError('run event must be a dictionary')
        event=row.get('event',{})
        checkpoint=(self.durability_profile=='PER_EVENT_FSYNC_V1' or row.get('kind')=='RUN_TERMINATED' or
            (row.get('kind')=='NEWTON_EVENT' and isinstance(event,dict) and
             event.get('type') in ('evaluation','attempt_start','attempt_terminal')))
        try:
            payload=self._serialize(row)+'\n'
        except ValueError:
            # Failed numerical trials can carry infinity as a diagnostic.
            # Preserve it as a typed JSON value; successful states and material
            # files continue to require finite numbers.
            if (row.get('kind')!='NEWTON_EVENT' or not isinstance(event,dict) or
                    event.get('outcome')=='accepted' or event.get('reason')=='converged'):
                raise
            from .faster_material import _json_value as encode_diagnostic
            tagged=encode_diagnostic(row)
            tagged['nonfinite_encoding']='IEEE754_TAGGED_V1'
            payload=self._serialize(tagged)+'\n'
        # Persist new material directory entries before writing the checkpoint
        # event that can acknowledge a completed or accepted evaluation.
        if checkpoint and self._pending_directory:self._directory_sync()
        start=time.monotonic_ns();self.events.write(payload);self.events.flush()
        self.metrics['write_flush_ns']+=time.monotonic_ns()-start
        if checkpoint:self._fsync(self.events.fileno(),'event')

    def close(self):
        if self._closed:return
        try:
            if self._pending_directory:self._directory_sync()
            self.events.flush();self._fsync(self.events.fileno(),'event')
        finally:
            self.events.close();self._closed=True
