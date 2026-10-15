"""Completed durable journals must equal the returned native stream."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from test_spin2014_authority_numerical_checker import checker


class AuthorityJournalChecks(unittest.TestCase):
    def test_completed_stream_and_truncated_or_mismatched_records(self):
        # Break caught: accepting missing, reordered, altered or wrong-input records.
        api = checker()
        identity = dict(request_identity="a"*64, context_identity="b"*64,
                        objective_identity="c"*64)
        events = [dict(type=kind, sequence=i, **identity)
                  for i,kind in enumerate(("attempt_start","attempt_terminal","run_terminal"))]
        run = dict(events=events, request_identity=identity["request_identity"],
                   context_identity=identity["context_identity"])
        expected = {"expected_"+key:value for key,value in identity.items()}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"journal.jsonl"
            def write(rows):
                path.write_text("".join(json.dumps(row)+"\n" for row in rows))
            rows = [dict(schema="spin2014-solver-event-journal-v1", sequence=i,event=event)
                    for i,event in enumerate(events)]
            write(rows)
            result = api.check_authority_event_journal(path, run=run, **expected)
            self.assertTrue(result["completed_stream_matches"])
            self.assertFalse(result["sealed_artifact"])
            variants = [rows[:-1], rows+[rows[-1]], list(reversed(rows))]
            changed = copy.deepcopy(rows);changed[0]["event"]["context_identity"]="d"*64
            variants.append(changed)
            for variant in variants:
                write(variant)
                with self.assertRaises(ValueError):
                    api.check_authority_event_journal(path, run=run, **expected)
            write(rows)
            path.write_text(path.read_text().rstrip("\n"))
            with self.assertRaises(ValueError):
                api.check_authority_event_journal(path, run=run, **expected)


if __name__ == "__main__":
    unittest.main()
