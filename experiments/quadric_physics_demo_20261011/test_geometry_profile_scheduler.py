import importlib.util
from pathlib import Path
import subprocess
import sys
import time
import unittest


class SchedulerTests(unittest.TestCase):
    def module(self):
        self.assertIsNotNone(importlib.util.find_spec('launch_native_geometry_profile'),
                             'serial diagnostic scheduler is missing')
        import launch_native_geometry_profile
        return launch_native_geometry_profile

    def test_failure_resumes_only_dispatcher_owned_by_this_context(self):
        module = self.module()
        child = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(30)', 'profile-test-dispatcher'])
        try:
            with self.assertRaisesRegex(RuntimeError, 'deliberate'):
                with module.paused_dispatcher(child.pid, 'profile-test-dispatcher'):
                    self.assertEqual(module.process_state(child.pid), 'T')
                    raise RuntimeError('deliberate')
            self.assertNotEqual(module.process_state(child.pid), 'T')
        finally:
            child.terminate(); child.wait()

    def test_wrong_identity_cannot_stop_another_process(self):
        module = self.module()
        child = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(30)', 'profile-test-dispatcher'])
        try:
            with self.assertRaisesRegex(RuntimeError, 'IDENTITY'):
                with module.paused_dispatcher(child.pid, 'not-this-process'):
                    self.fail('wrong process was accepted')
            self.assertNotEqual(module.process_state(child.pid), 'T')
        finally:
            child.terminate(); child.wait()


if __name__ == '__main__':
    unittest.main()
