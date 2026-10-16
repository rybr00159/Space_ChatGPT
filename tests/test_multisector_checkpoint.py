"""Regression: SQLite backup archives must include the destination WAL state."""
from contextlib import closing
from pathlib import Path
import importlib.util,json,sqlite3,tempfile,unittest


ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('checkpoint_under_test',ROOT/'scripts/multisector_checkpoint.py')
CHECKPOINT=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(CHECKPOINT)


class CheckpointWALTest(unittest.TestCase):
    def test_repeated_backup_contains_new_rows_and_restores_without_sidecars(self):
        previous_cache,previous_out=CHECKPOINT.CACHE,CHECKPOINT.OUT
        with tempfile.TemporaryDirectory(prefix='tess_checkpoint_wal_')as temporary:
            folder=Path(temporary);CHECKPOINT.CACHE=folder/'cache';CHECKPOINT.OUT=folder/'tracked'
            CHECKPOINT.CACHE.mkdir()
            try:
                # Deliberately keep a source writer open with active WAL pages.
                with closing(sqlite3.connect(CHECKPOINT.CACHE/'survey.sqlite'))as writer:
                    writer.execute('PRAGMA journal_mode=WAL')
                    writer.execute('CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT)')
                    writer.execute("INSERT INTO metadata VALUES('methods','{}')")
                    writer.execute('CREATE TABLE observations(tic INTEGER,sector INTEGER,status TEXT,result_json TEXT)')
                    writer.execute("INSERT INTO observations VALUES(1,83,'synthetic_backup_test',NULL)")
                    writer.commit();CHECKPOINT.freeze()
                    CHECKPOINT.restore(folder/'first.sqlite')
                    writer.executemany("INSERT INTO observations VALUES(?,86,'synthetic_backup_test',NULL)",[(2,),(3,)])
                    writer.commit();CHECKPOINT.freeze()
                    restored=folder/'second.sqlite';CHECKPOINT.restore(restored)
                    with closing(sqlite3.connect(restored))as reader:
                        self.assertEqual(reader.execute('SELECT tic FROM observations ORDER BY tic').fetchall(),[(1,),(2,),(3,)])
                        self.assertEqual(reader.execute('PRAGMA journal_mode').fetchone()[0],'delete')
                    self.assertFalse(restored.with_name(restored.name+'-wal').exists())
                    self.assertEqual(json.loads((CHECKPOINT.OUT/'manifest.json').read_text())['completed_observations'],3)
            finally:CHECKPOINT.CACHE,CHECKPOINT.OUT=previous_cache,previous_out


if __name__=='__main__':unittest.main()
