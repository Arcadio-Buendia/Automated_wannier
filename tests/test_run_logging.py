import contextlib
import io
import json
from datetime import datetime
from pathlib import Path
import os
import time
from unittest.mock import patch
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"workflow"))
from extensions.run_logging import log, timestamp
from wannier_workflow import Workflow


class LoggingTests(unittest.TestCase):
    def test_timezone_and_progress(self):
        self.assertEqual(datetime.fromisoformat(timestamp()).utcoffset(), datetime.now().astimezone().utcoffset())
        output = io.StringIO()
        log('hello', file=output)
        prefix, message = output.getvalue().split('] ', 1)
        self.assertEqual(datetime.fromisoformat(prefix[1:]).utcoffset(), datetime.now().astimezone().utcoffset())
        self.assertEqual(message, 'hello\n')

    @unittest.skipUnless(hasattr(time, 'tzset'), 'Requires POSIX timezone support')
    def test_machine_timezone_is_used(self):
        try:
            with patch.dict(os.environ, {'TZ': 'UTC-5:30'}):
                time.tzset()
                self.assertTrue(timestamp().endswith('+05:30'))
        finally:
            time.tzset()

    def test_stage_events_preserve_raw_output_and_record_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            runner = Workflow.__new__(Workflow)
            runner.work = root
            runner.config = {'commands': {'pw': [sys.executable, '-c', 'print("JOB DONE.")']}}
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(runner.run_command('pw', root), 0)
            self.assertEqual((root/'pw.out').read_text(), 'JOB DONE.\n')
            runner.config['commands']['pw'] = [str(root/'missing_executable')]
            with contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaises(FileNotFoundError):
                    runner.run_command('pw', root)
            events = [json.loads(line) for line in (root/'stage_events.jsonl').read_text().splitlines()]
            self.assertEqual(len(events), 4)
            self.assertNotIn('finished_at', events[0])
            self.assertEqual(events[1]['returncode'], 0)
            self.assertGreaterEqual(events[1]['elapsed_seconds'], 0)
            self.assertIn('FileNotFoundError', events[3]['error'])
            for event in events:
                self.assertEqual(datetime.fromisoformat(event['started_at']).utcoffset(), datetime.now().astimezone().utcoffset())
