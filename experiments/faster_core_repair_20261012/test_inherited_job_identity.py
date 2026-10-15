import hashlib
import pytest


def test_namespaced_dependencies_preserve_reserved_entries_and_reject_tampering(tmp_path):
    from inherited_job_identity import inherited_files
    source = tmp_path / 'source'
    source.write_text('bound')
    item = {'path': str(source), 'sha256': hashlib.sha256(b'bound').hexdigest()}
    result = inherited_files({'__job_worker': item, 'worker': item})
    assert len(result) == 2
    assert all(not key.startswith('__job_') for key in result)
    assert list(result.values()) == [source, source]
    source.write_text('changed')
    with pytest.raises(ValueError, match='DEPENDENCY_CHANGED'):
        inherited_files({'__job_worker': item})
