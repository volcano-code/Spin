# Task5a controller verification at f2998d0

Main work: explicit paper Givens formulas and actual joint fill/angle integration. Independent review approved bounded task, two deferred Minor notes, no Critical/Important.

## Fresh observed checks

- Full configured CTest189/189 actualpasses, exit0,260.62seconds, maximumRSS4253708KiB. NoNotRun/Skipped/Disabled. FullrawlogSHA d30e2b58ffb5516df573de264871b098393ed554143f7f6541f3516976abb66b; consoleSHA a5671140457ee4147714ecebceeceec3511fd0558ac7e348242902ff5d403d14. Logs preserved beside task-5a-report.md.
- Controller postcommit CTest4/4(givens,slqp,end-to-end,python), actual8Pythonmethods,0.12seconds. Raw task-5a-controller-postcommit.log retained.
- `sha256sum -c build/strict-v2-spin/task5a-final-sha256-20260908.log`: all23entriesOK. CurrentDSO4052d643674d5538851f4aa2f711b2b960ecdf1a26c5b6242265aa80aea0317a. ActualCMakeCache2014/Faster/exact/PythonON, viewer/PyBulletOFF.
- Synthetic joint fixture: initialphi1.7, fouracceptedsteps, finalphi1.5707963267947704, equality∞1.4213146065523828e-11. Maximumfourthorderderivativeerror1.2127399778849577e-10 withfixedthreeh and originalsecondorderdiagnosticsretained. Adaptiveanglecarryobservedtwo passes. This is not a paper-model result or physical simulation.
- MaintrackedSHA93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e unchanged. Stricttrackedtreeclean atimplementationhandoff. Source/testdiffwhitespacecheckclean; preservedrawconsolewhitespaceisnotrewritten.

## Existing result freshness and narrow evidence refresh

The previous Task4c3a candidate bundle is preserved at build/task-4c3a-candidate-ae1OAozK/material. Actualcurrentloader rejects it with `stale implementation observation` after nativeDSO changes. That is correct provenance rejection, not evidence of a newly invalid geometry. No attempt to silently reseal/relabel it was made. CurrentmanifestSHA2e314f2766f24514c07a03f99a2080a246c1503a8c72b9cda915fdfafb002f82 recorded only as observation, not an unmeasured before/after comparison.

Task7arun3/4 are preserved and nowSTALE_RESULT because declared `cmake_configuration` input(CMakeLists.txt) changed with testregistration. After review/codefreeze, controller reran the unchanged exact Section6.6 runner twice in NEW directories:

- build/task-5a-refresh-faster-section-6-6-run-1
- build/task-5a-refresh-faster-section-6-6-run-2

Command for each N:

```sh
env PYTHONPATH=python /home/lin/anaconda3/bin/python3 experiments/spin_it_reproduction/strict/run_faster_section_6_6.py --spec experiments/spin_it_reproduction/sources/faster_section_6_6_spec.json --paper build/strict-paper-sources/sif-final.pdf --output-directory build/task-5a-refresh-faster-section-6-6-run-N
```

BothCURRENT. Controller separately read all5input+2output files for eachrun, checkedactualsizes/SHAagainstrecords, and comparedresult/reportbytes acrossruns. ScientificSHA4a8e15c902415d85e381c9290f56398cd4f8fbd89d0753360d708e380d3a4754; resultSHA4f6e763cb1b681babe4dfddfb2dcae1743a7493897a70234a6ba9e1ce2ab72fb; reportSHA4a1cd2efc0febb98974b5d70a2a6e9961f3dea2527f6587b77c9a503856f5d98. Valuesunchangedfromhistoricalproofs; provenancecurrentinputdigest1e72c66a3b40a1a86bcf089b426adb55ea903080eced3f8490349c68f5791878. CurrentCMakeListsSHA6b9234409fb177266d828504ab82518298903150c46143bc5d956d9dc63db4c0.

Originalspec/PDFhashes still053f896c2b10ece89b136bfba24cd480bf44cd0da99c0661e6e9bd34495c15fa / e7fc6395362ef7e63a703851176afd81d8e953bae4e8080279d7e244b3c51fcc. No oldscientificfile changed. These are the same two published1Ddegenerateexamples, not two newcases or3D/CAD/physical results. Their purePythonproofdoesnotattestnativebuild. Alloriginal3Dcase/fullmethodreleasegates remainopen.
