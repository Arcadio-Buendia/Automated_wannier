import contextlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"workflow"))
from extensions.run_logging import log, timestamp
from wannier_workflow import Workflow


class LoggingTests(unittest.TestCase):
    def test_timezone_and_progress(self):
        self.assertEqual(datetime.fromisoformat(timestamp()).tzinfo, timezone.utc)
        output = io.StringIO()
        log('hello', file=output)
        prefix, message = output.getvalue().split('] ', 1)
        self.assertEqual(datetime.fromisoformat(prefix[1:]).tzinfo, timezone.utc)
        self.assertEqual(message, 'hello\n')

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
                self.assertEqual(datetime.fromisoformat(event['started_at']).tzinfo, timezone.utc)
