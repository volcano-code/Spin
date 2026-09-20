# Static archive member-byte audit

This is a deliberately narrow, historical observation. It binds the member
payloads of two ordinary GNU `ar` archives to the object bytes in independent,
successful compiler receipts. It is not a build execution, source, linker,
runtime, solver, paper, or physics attestation.

## Historical inputs

The preserved input build is
`build/strict-v2-observed-zf5mMW`, recorded against source revision
`7057c8120997c7e0d70b8d1a11833822be3c9ed9`, not this dirty worktree. Its
source-content digest was
`afd3ee994984618c7df9a5047466984ccc7a61b207e1a6ec712979815cbd8355` before
and after its recorded postchecks.

The verifier independently read and bound these input bytes:

- `build-inventory-v3.json`:
  `13105c0a871cbb2cdfa7951435122c1eb47a09cc94688ca3c4cd7a46b7c8fa3f`
- `ninja-link-evidence-v2.json`:
  `16c7608ed9bbfbee9e7671d17c5202c88df5f0e239e3b8993f4ac2c2fab6f875`
- generated `ninja-commands.stdout`:
  `d641be323253d09b32f88650f90f83b2905ed7cdf084f589b73afbd0fc2c19f3`

`/usr/bin/ar` was itself hashed as
`534681ac11c18868cfc4fdf98770aa0ba8973eedc90c231e94e6ba96e1a04f27`
(55,792 bytes; resolved executable `/usr/bin/x86_64-linux-gnu-ar`).

## Result

Both target archive results are `OBSERVED_ARCHIVE_MEMBER_BYTES_MATCH`:

| Archive | Archive SHA-256 | Archive bytes | Members | Bound member bytes | Result JSON SHA-256 |
| --- | --- | ---: | ---: | ---: | --- |
| `build/lib/libhighs.a` | `4e90d129307e57b0045ae1b057b40bc2f3a7c6b358aafbb9579f78d4309c7b77` | 283,230,962 | 181 | 282,859,240 | `efb19db8314e057ba0b50a9695abfc2ddbc260be1fd7e84e617b650e95cd63c6` |
| `build/libmake_it_stand_core.a` | `4091671a2057b985bbba72ea357e67960c3c1ede7c3a5cd0bd4abc2bfe4c6e13` | 561,039,820 | 57 | 560,056,944 | `94e7be99811333749c88d640b8ed5d2fedfd9327ead1dcccd0944e5f415eea40` |

Machine-readable outputs are intentionally ignored under
`build/strict-v2-archive-member-audit-20260921-line-bound/`; the historical build and
all of its prior evidence were only read.

## Rerun

From the repository root, use a new output directory (the CLI refuses to
overwrite a result):

```sh
PYTHONPATH=python python3 -m make_it_stand.verification.source_archive_member_audit \
  --inventory build/strict-v2-observed-zf5mMW/build-inventory-v3.json \
  --ninja-evidence build/strict-v2-observed-zf5mMW/ninja-link-evidence-v2.json \
  --archive build/strict-v2-observed-zf5mMW/build/lib/libhighs.a \
  --output build/new-archive-audit/libhighs.member-audit.json
```

Repeat with `--archive build/strict-v2-observed-zf5mMW/build/libmake_it_stand_core.a`.
The verifier admits only the captured GNU/CMake `ar qc … && ranlib` grammar,
ordinary `!<arch>\n` archives, unique safe member basenames, and one matching
`target_compilation`/`OBSERVED_SUCCESS` compiler receipt per object argument.
It hashes raw receipts, current objects, member payloads, and the archive both
before and after inspection. Each static-archive structured record must also
match its declared 1-based, freshly rehashed raw `ninja-commands.stdout` line;
invalid, duplicate, or detached line records fail closed. Thin, malformed,
ambiguous, duplicate, stale, and mutating inputs fail closed.

## Explicit limitations

The generated Ninja command list is **not** a trace proving that an archive
command executed. This work does not close implicit linker inputs, assembler or
compiler plugins, external prebuilt libraries, source/build provenance, or
loaded runtime bytes. It never issues `VERIFIED`, solver acceptance, paper
acceptance, or physics acceptance.
