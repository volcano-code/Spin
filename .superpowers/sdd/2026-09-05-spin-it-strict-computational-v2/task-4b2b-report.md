# Task 4b2b report — static archive member binding

## Scope and historical input

Implemented a bounded verifier at
`python/make_it_stand/verification/source_archive_member_audit.py`, with real
GNU-`ar` behavioral fixtures and a CTest registration. The audited input is the
preserved Task 4b2a build `build/strict-v2-observed-zf5mMW`, recorded against
source revision `7057c8120997c7e0d70b8d1a11833822be3c9ed9`; it is not a build
of the current dirty source tree. No file under that historical evidence
directory was modified or regenerated.

Input SHA-256 values: inventory
`13105c0a871cbb2cdfa7951435122c1eb47a09cc94688ca3c4cd7a46b7c8fa3f`,
Ninja evidence
`16c7608ed9bbfbee9e7671d17c5202c88df5f0e239e3b8993f4ac2c2fab6f875`, and
raw generated-command output
`d641be323253d09b32f88650f90f83b2905ed7cdf084f589b73afbd0fc2c19f3`.
The historical source content digest was
`afd3ee994984618c7df9a5047466984ccc7a61b207e1a6ec712979815cbd8355` before
and after recorded postchecks.

## TDD evidence

RED command:

```sh
PYTHONPATH=python pytest -q tests/python/test_source_archive_member_audit.py
```

Initial output: `11 failed`: ten intended behavioral cases failed with
`ModuleNotFoundError: No module named
'make_it_stand.verification.source_archive_member_audit'`, while the thin
fixture initially tried to convert an existing ordinary archive and failed in
GNU `ar` setup. That fixture setup was corrected before implementation; no
production verifier existed during the RED run.

GREEN command:

```sh
PYTHONPATH=python PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 pytest -q tests/python/test_source_archive_member_audit.py
```

Output: `12 passed in 0.28s`. The real fixtures cover exact member-byte
matching, altered same-name payloads, thin archives, duplicate basenames,
missing/extra/reordered/unsafe members, malformed receipt and object identity,
digest mismatch, and mutation after member inspection.

Registered CTest command:

```sh
cmake -S . -B build/strict-v2-spin
ctest --test-dir build/strict-v2-spin -R '^source-archive-member-audit-contract$' --output-on-failure
```

Output: `1/1 Test #139: source-archive-member-audit-contract ... Passed`
and `100% tests passed, 0 tests failed out of 1`.

## Historical archive results

The following read-only commands wrote new machine-readable outputs below the
ignored task-specific directory `build/strict-v2-archive-member-audit-20260921-final/`:

```sh
PYTHONPATH=python python3 -m make_it_stand.verification.source_archive_member_audit \
  --inventory build/strict-v2-observed-zf5mMW/build-inventory-v3.json \
  --ninja-evidence build/strict-v2-observed-zf5mMW/ninja-link-evidence-v2.json \
  --archive build/strict-v2-observed-zf5mMW/build/lib/libhighs.a \
  --output build/strict-v2-archive-member-audit-20260921-final/libhighs.member-audit.json
PYTHONPATH=python python3 -m make_it_stand.verification.source_archive_member_audit \
  --inventory build/strict-v2-observed-zf5mMW/build-inventory-v3.json \
  --ninja-evidence build/strict-v2-observed-zf5mMW/ninja-link-evidence-v2.json \
  --archive build/strict-v2-observed-zf5mMW/build/libmake_it_stand_core.a \
  --output build/strict-v2-archive-member-audit-20260921-final/libmake_it_stand_core.member-audit.json
```

Both exit statuses were zero and each result status is
`OBSERVED_ARCHIVE_MEMBER_BYTES_MATCH`.

| Archive | Archive SHA-256 / bytes | Members / member bytes | Result SHA-256 |
| --- | --- | --- | --- |
| `libhighs.a` | `4e90d129307e57b0045ae1b057b40bc2f3a7c6b358aafbb9579f78d4309c7b77` / 283,230,962 | 181 / 282,859,240 | `efb19db8314e057ba0b50a9695abfc2ddbc260be1fd7e84e617b650e95cd63c6` |
| `libmake_it_stand_core.a` | `4091671a2057b985bbba72ea357e67960c3c1ede7c3a5cd0bd4abc2bfe4c6e13` / 561,039,820 | 57 / 560,056,944 | `94e7be99811333749c88d640b8ed5d2fedfd9327ead1dcccd0944e5f415eea40` |

The verifier recorded `/usr/bin/ar` as resolved
`/usr/bin/x86_64-linux-gnu-ar`, 55,792 bytes, SHA-256
`534681ac11c18868cfc4fdf98770aa0ba8973eedc90c231e94e6ba96e1a04f27`.

## Files changed

- `python/make_it_stand/verification/source_archive_member_audit.py`
- `tests/python/test_source_archive_member_audit.py`
- `CMakeLists.txt` (one additive CTest registration hunk only)
- `docs/source_archive_member_audit.md`
- this report

## Limitations

Ninja command evidence is generated-command evidence, not an archive-command
execution trace. The result binds only archive member payloads to independently
rehashable compiler object receipts. It does not grant build-source, compiler or
assembler plugin, linker implicit-input, runtime-memory, solver, paper, or
physics acceptance; all of those remain `NOT_VERIFIED`.
