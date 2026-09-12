cmake_minimum_required(VERSION 3.16)

if(NOT DEFINED TEST_BINARY_ROOT OR TEST_BINARY_ROOT STREQUAL "")
    message(FATAL_ERROR "TEST_BINARY_ROOT is required")
endif()
if(NOT IS_ABSOLUTE "${TEST_BINARY_ROOT}" OR NOT IS_DIRECTORY "${TEST_BINARY_ROOT}")
    message(FATAL_ERROR "TEST_BINARY_ROOT must be an existing absolute directory")
endif()
if(NOT DEFINED DEPENDENCY_SCRIPT OR NOT EXISTS "${DEPENDENCY_SCRIPT}")
    message(FATAL_ERROR "DEPENDENCY_SCRIPT must name the dependency script under test")
endif()

string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef fixture_suffix)
set(FIXTURE_ROOT "${TEST_BINARY_ROOT}/spin-2014-isolation-fixture-${fixture_suffix}")
set(SUCCESS_SOURCE "${FIXTURE_ROOT}/success-source")
set(SUCCESS_ROOT "${FIXTURE_ROOT}/success-invocations")
file(MAKE_DIRECTORY
    "${SUCCESS_SOURCE}/tests/cmake"
    "${SUCCESS_SOURCE}/src/spin2014"
    "${SUCCESS_ROOT}"
)
configure_file(
    "${DEPENDENCY_SCRIPT}"
    "${SUCCESS_SOURCE}/tests/cmake/assert_spin_it_2014_dependency.cmake"
    COPYONLY
)
if(DEFINED ALLOCATOR_HELPER AND EXISTS "${ALLOCATOR_HELPER}")
    configure_file(
        "${ALLOCATOR_HELPER}"
        "${SUCCESS_SOURCE}/tests/cmake/allocate_retained_test_directory.cmake"
        COPYONLY
    )
endif()
file(WRITE "${SUCCESS_SOURCE}/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.16)
project(spin_2014_dependency_fixture LANGUAGES CXX)
option(MAKE_IT_STAND_BUILD_SPIN_IT_2014 "" OFF)
option(LIBIGL_SKIP_DOWNLOAD "" OFF)
if(NOT LIBIGL_SKIP_DOWNLOAD)
    message(FATAL_ERROR "fixture dependency downloads were not disabled")
endif()
if(MAKE_IT_STAND_BUILD_SPIN_IT_2014)
    add_executable(
        spin-it-2014-audited
        src/spin2014/AdaptiveVoxelSpinSolver.cpp
    )
endif()
]=])
file(WRITE "${SUCCESS_SOURCE}/src/spin2014/AdaptiveVoxelSpinSolver.cpp" [=[
#include <iostream>
int main() {
    std::cout << "fixture audited pass";
    return 0;
}
]=])

foreach(old_name IN ITEMS spin-2014-disabled-contract spin-2014-enabled-contract)
    file(MAKE_DIRECTORY "${SUCCESS_SOURCE}/build/${old_name}")
    file(WRITE "${SUCCESS_SOURCE}/build/${old_name}/sentinel.txt" "preserve-${old_name}\n")
endforeach()

function(run_success_contract output_var)
    execute_process(
        COMMAND
            "${CMAKE_COMMAND}"
            "-DMAKE_IT_STAND_DEPENDENCY_TEST_BINARY_ROOT=${SUCCESS_ROOT}"
            -P "${SUCCESS_SOURCE}/tests/cmake/assert_spin_it_2014_dependency.cmake"
        RESULT_VARIABLE contract_status
        OUTPUT_VARIABLE contract_stdout
        ERROR_VARIABLE contract_stderr
    )
    if(NOT contract_status EQUAL 0)
        message(FATAL_ERROR
            "fixture dependency contract failed (${contract_status})\n"
            "stdout:\n${contract_stdout}\nstderr:\n${contract_stderr}"
        )
    endif()
    set(contract_output "${contract_stdout}\n${contract_stderr}")
    string(REGEX MATCH
        "SPIN_IT_2014_DEPENDENCY_RETAINED_DIR=([^\n\r]+)"
        retained_match
        "${contract_output}"
    )
    if(retained_match)
        set(${output_var} "${CMAKE_MATCH_1}" PARENT_SCOPE)
    else()
        set(${output_var} "" PARENT_SCOPE)
    endif()
endfunction()

run_success_contract(first_run)
foreach(old_name IN ITEMS spin-2014-disabled-contract spin-2014-enabled-contract)
    set(sentinel "${SUCCESS_SOURCE}/build/${old_name}/sentinel.txt")
    if(NOT EXISTS "${sentinel}")
        message(FATAL_ERROR "dependency contract deleted historical sentinel: ${sentinel}")
    endif()
    file(READ "${sentinel}" sentinel_text)
    if(NOT sentinel_text STREQUAL "preserve-${old_name}\n")
        message(FATAL_ERROR "dependency contract altered historical sentinel: ${sentinel}")
    endif()
endforeach()
if(first_run STREQUAL "")
    message(FATAL_ERROR "contract did not print its retained directory before work")
endif()
if(NOT IS_DIRECTORY "${first_run}")
    message(FATAL_ERROR "printed retained directory does not exist: ${first_run}")
endif()
file(GLOB_RECURSE first_files RELATIVE "${first_run}" "${first_run}/*")
list(SORT first_files)
set(first_snapshot "")
foreach(relative_file IN LISTS first_files)
    if(NOT IS_DIRECTORY "${first_run}/${relative_file}")
        file(SHA256 "${first_run}/${relative_file}" file_hash)
        string(APPEND first_snapshot "${relative_file}=${file_hash}\n")
    endif()
endforeach()

run_success_contract(second_run)
if(first_run STREQUAL second_run)
    message(FATAL_ERROR "two dependency invocations reused ${first_run}")
endif()
file(GLOB_RECURSE first_files_after RELATIVE "${first_run}" "${first_run}/*")
list(SORT first_files_after)
set(first_snapshot_after "")
foreach(relative_file IN LISTS first_files_after)
    if(NOT IS_DIRECTORY "${first_run}/${relative_file}")
        file(SHA256 "${first_run}/${relative_file}" file_hash)
        string(APPEND first_snapshot_after "${relative_file}=${file_hash}\n")
    endif()
endforeach()
if(NOT first_snapshot STREQUAL first_snapshot_after)
    message(FATAL_ERROR "second dependency invocation altered the first retained run")
endif()

function(expect_rejected_root label root_argument)
    execute_process(
        COMMAND
            "${CMAKE_COMMAND}"
            "-DMAKE_IT_STAND_DEPENDENCY_TEST_BINARY_ROOT=${root_argument}"
            -P "${SUCCESS_SOURCE}/tests/cmake/assert_spin_it_2014_dependency.cmake"
        RESULT_VARIABLE rejected_status
        OUTPUT_VARIABLE rejected_stdout
        ERROR_VARIABLE rejected_stderr
    )
    if(rejected_status EQUAL 0)
        message(FATAL_ERROR "${label} dependency root was accepted")
    endif()
endfunction()

execute_process(
    COMMAND "${CMAKE_COMMAND}" -P "${SUCCESS_SOURCE}/tests/cmake/assert_spin_it_2014_dependency.cmake"
    RESULT_VARIABLE missing_status
    OUTPUT_VARIABLE missing_stdout
    ERROR_VARIABLE missing_stderr
)
if(missing_status EQUAL 0)
    message(FATAL_ERROR "missing dependency root was accepted")
endif()
expect_rejected_root("relative" "relative-build-root")
expect_rejected_root("nonexistent" "${FIXTURE_ROOT}/does-not-exist")
expect_rejected_root("source" "${SUCCESS_SOURCE}")
if(EXISTS "${FIXTURE_ROOT}/does-not-exist")
    message(FATAL_ERROR "rejected nonexistent root was created")
endif()
if(EXISTS "${SUCCESS_SOURCE}/.spin-2014-dependency-contract-allocator.lock"
   OR EXISTS "${SUCCESS_SOURCE}/spin-2014-dependency-contract-runs")
    message(FATAL_ERROR "rejected source root was mutated")
endif()
foreach(old_name IN ITEMS spin-2014-disabled-contract spin-2014-enabled-contract)
    set(sentinel "${SUCCESS_SOURCE}/build/${old_name}/sentinel.txt")
    file(READ "${sentinel}" sentinel_text)
    if(NOT sentinel_text STREQUAL "preserve-${old_name}\n")
        message(FATAL_ERROR "rejected root invocation altered historical sentinel: ${sentinel}")
    endif()
endforeach()

set(ESCAPE_ROOT "${FIXTURE_ROOT}/escape-root")
set(ESCAPE_OUTSIDE "${FIXTURE_ROOT}/escape-outside")
file(MAKE_DIRECTORY "${ESCAPE_ROOT}" "${ESCAPE_OUTSIDE}")
file(CREATE_LINK
    "${ESCAPE_OUTSIDE}"
    "${ESCAPE_ROOT}/spin-2014-dependency-contract-runs"
    SYMBOLIC
    RESULT escape_link_status
)
if(escape_link_status)
    message(FATAL_ERROR "could not create symlink escape fixture: ${escape_link_status}")
endif()
expect_rejected_root("symlink-escaping" "${ESCAPE_ROOT}")
file(GLOB escape_outside_files "${ESCAPE_OUTSIDE}/*")
if(escape_outside_files)
    message(FATAL_ERROR "symlink-escaping root mutated outside destination")
endif()

if(NOT DEFINED ALLOCATOR_HELPER OR NOT EXISTS "${ALLOCATOR_HELPER}")
    message(FATAL_ERROR "allocator helper is missing")
endif()
set(CONCURRENT_SOURCE "${FIXTURE_ROOT}/concurrent-source")
set(CONCURRENT_BUILD "${FIXTURE_ROOT}/concurrent-build")
set(CONCURRENT_ROOT "${FIXTURE_ROOT}/concurrent-root")
file(MAKE_DIRECTORY "${CONCURRENT_SOURCE}" "${CONCURRENT_ROOT}")
file(WRITE "${CONCURRENT_SOURCE}/allocate.cmake" [=[
include("${ALLOCATOR_HELPER}")
make_it_stand_allocate_retained_test_directory("${ALLOCATION_ROOT}" allocated_dir)
file(WRITE "${RESULT_FILE}" "${allocated_dir}\n")
]=])
file(WRITE "${CONCURRENT_SOURCE}/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.16)
project(concurrent_allocator_fixture NONE)
foreach(worker IN ITEMS a b)
    add_custom_target("allocate-${worker}"
        COMMAND "${CMAKE_COMMAND}"
            "-DALLOCATOR_HELPER=${ALLOCATOR_HELPER}"
            "-DALLOCATION_ROOT=${ALLOCATION_ROOT}"
            "-DRESULT_FILE=${ALLOCATION_ROOT}/concurrent-${worker}.txt"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/allocate.cmake"
    )
endforeach()
add_custom_target(allocate-all DEPENDS allocate-a allocate-b)
]=])
execute_process(
    COMMAND "${CMAKE_COMMAND}"
        -S "${CONCURRENT_SOURCE}"
        -B "${CONCURRENT_BUILD}"
        -G Ninja
        "-DALLOCATOR_HELPER=${ALLOCATOR_HELPER}"
        "-DALLOCATION_ROOT=${CONCURRENT_ROOT}"
    RESULT_VARIABLE concurrent_configure_status
    OUTPUT_VARIABLE concurrent_configure_stdout
    ERROR_VARIABLE concurrent_configure_stderr
)
if(NOT concurrent_configure_status EQUAL 0)
    message(FATAL_ERROR "concurrent fixture configure failed: ${concurrent_configure_stderr}")
endif()
execute_process(
    COMMAND "${CMAKE_COMMAND}" --build "${CONCURRENT_BUILD}" --target allocate-all -j2
    RESULT_VARIABLE concurrent_status
    OUTPUT_VARIABLE concurrent_stdout
    ERROR_VARIABLE concurrent_stderr
)
if(NOT concurrent_status EQUAL 0)
    message(FATAL_ERROR
        "concurrent allocation failed (${concurrent_status})\n"
        "stdout:\n${concurrent_stdout}\nstderr:\n${concurrent_stderr}"
    )
endif()
file(READ "${CONCURRENT_ROOT}/concurrent-a.txt" concurrent_a)
file(READ "${CONCURRENT_ROOT}/concurrent-b.txt" concurrent_b)
string(STRIP "${concurrent_a}" concurrent_a)
string(STRIP "${concurrent_b}" concurrent_b)
if(concurrent_a STREQUAL concurrent_b)
    message(FATAL_ERROR "concurrent allocations shared ${concurrent_a}")
endif()
if(NOT IS_DIRECTORY "${concurrent_a}" OR NOT IS_DIRECTORY "${concurrent_b}")
    message(FATAL_ERROR "concurrent allocation result is missing")
endif()

set(LIVE_SOURCE "${FIXTURE_ROOT}/live-source")
set(LIVE_ROOT "${FIXTURE_ROOT}/live-invocation")
set(LIVE_DRIVER_SOURCE "${FIXTURE_ROOT}/live-driver-source")
set(LIVE_DRIVER_BUILD "${FIXTURE_ROOT}/live-driver-build")
file(COPY "${SUCCESS_SOURCE}/" DESTINATION "${LIVE_SOURCE}")
file(MAKE_DIRECTORY "${LIVE_ROOT}" "${LIVE_DRIVER_SOURCE}")
file(WRITE "${LIVE_SOURCE}/src/spin2014/AdaptiveVoxelSpinSolver.cpp" [=[
#include <chrono>
#include <iostream>
#include <thread>
int main() {
    std::cout << "live stdout before sleep" << std::endl;
    std::cerr << "live stderr before sleep" << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(5));
    return 0;
}
]=])
file(WRITE "${LIVE_DRIVER_SOURCE}/observe.cmake" [=[
set(observed FALSE)
foreach(attempt RANGE 1 40)
    file(GLOB live_runs "${LIVE_ROOT}/spin-2014-dependency-contract-runs/run-*")
    foreach(live_run IN LISTS live_runs)
        if(EXISTS "${live_run}/audited.status"
           AND EXISTS "${live_run}/audited.stdout.log"
           AND EXISTS "${live_run}/audited.stderr.log")
            file(READ "${live_run}/audited.status" live_status)
            file(READ "${live_run}/audited.stdout.log" live_stdout)
            file(READ "${live_run}/audited.stderr.log" live_stderr)
            string(STRIP "${live_status}" live_status)
            if(live_status STREQUAL "RUNNING"
               AND live_stdout MATCHES "live stdout before sleep"
               AND live_stderr MATCHES "live stderr before sleep")
                file(WRITE "${LIVE_ROOT}/observed-live.txt" "${live_run}\n")
                set(observed TRUE)
                break()
            endif()
        endif()
    endforeach()
    if(observed)
        break()
    endif()
    execute_process(COMMAND "${CMAKE_COMMAND}" -E sleep 0.25)
endforeach()
if(NOT observed)
    message(FATAL_ERROR "child logs and RUNNING status were not visible before child exit")
endif()
]=])
file(WRITE "${LIVE_DRIVER_SOURCE}/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.16)
project(live_log_fixture NONE)
add_custom_target(run-contract
    COMMAND "${CMAKE_COMMAND}"
        "-DMAKE_IT_STAND_DEPENDENCY_TEST_BINARY_ROOT=${LIVE_ROOT}"
        -P "${LIVE_SOURCE}/tests/cmake/assert_spin_it_2014_dependency.cmake"
)
add_custom_target(observe-contract
    COMMAND "${CMAKE_COMMAND}"
        "-DLIVE_ROOT=${LIVE_ROOT}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/observe.cmake"
)
add_custom_target(run-and-observe ALL DEPENDS run-contract observe-contract)
]=])
execute_process(
    COMMAND "${CMAKE_COMMAND}"
        -S "${LIVE_DRIVER_SOURCE}"
        -B "${LIVE_DRIVER_BUILD}"
        -G Ninja
        "-DLIVE_ROOT=${LIVE_ROOT}"
        "-DLIVE_SOURCE=${LIVE_SOURCE}"
    RESULT_VARIABLE live_configure_status
    OUTPUT_VARIABLE live_configure_stdout
    ERROR_VARIABLE live_configure_stderr
)
if(NOT live_configure_status EQUAL 0)
    message(FATAL_ERROR "live-log fixture configure failed: ${live_configure_stderr}")
endif()
execute_process(
    COMMAND "${CMAKE_COMMAND}" --build "${LIVE_DRIVER_BUILD}" -j2
    RESULT_VARIABLE live_status
    OUTPUT_VARIABLE live_stdout
    ERROR_VARIABLE live_stderr
)
if(NOT live_status EQUAL 0)
    message(FATAL_ERROR
        "live-log fixture failed (${live_status})\n"
        "stdout:\n${live_stdout}\nstderr:\n${live_stderr}"
    )
endif()
if(NOT EXISTS "${LIVE_ROOT}/observed-live.txt")
    message(FATAL_ERROR "live-log fixture did not record pre-exit observation")
endif()

set(FAILURE_SOURCE "${FIXTURE_ROOT}/failure-source")
set(FAILURE_ROOT "${FIXTURE_ROOT}/failure-invocations")
file(COPY "${SUCCESS_SOURCE}/" DESTINATION "${FAILURE_SOURCE}")
file(MAKE_DIRECTORY "${FAILURE_ROOT}")
file(WRITE "${FAILURE_SOURCE}/src/spin2014/AdaptiveVoxelSpinSolver.cpp" [=[
#include <iostream>
int main() {
    std::cout << "fixture audited failure";
    std::cerr << "fixture failure stderr";
    return 23;
}
]=])
execute_process(
    COMMAND
        "${CMAKE_COMMAND}"
        "-DMAKE_IT_STAND_DEPENDENCY_TEST_BINARY_ROOT=${FAILURE_ROOT}"
        -P "${FAILURE_SOURCE}/tests/cmake/assert_spin_it_2014_dependency.cmake"
    RESULT_VARIABLE failure_status
    OUTPUT_VARIABLE failure_stdout
    ERROR_VARIABLE failure_stderr
)
if(failure_status EQUAL 0)
    message(FATAL_ERROR "failing audited child was reported as passing")
endif()
set(failure_output "${failure_stdout}\n${failure_stderr}")
string(REGEX MATCH
    "SPIN_IT_2014_DEPENDENCY_RETAINED_DIR=([^\n\r]+)"
    failure_retained_match
    "${failure_output}"
)
if(NOT failure_retained_match)
    message(FATAL_ERROR "failing contract did not expose retained evidence")
endif()
set(failure_run "${CMAKE_MATCH_1}")
foreach(required_file IN ITEMS
    configure-disabled.stdout.log
    configure-disabled.stderr.log
    configure-disabled.status
    configure-enabled.stdout.log
    configure-enabled.stderr.log
    configure-enabled.status
    build-enabled.stdout.log
    build-enabled.stderr.log
    build-enabled.status
    audited.stdout.log
    audited.stderr.log
    audited.status
    contract.status
)
    if(NOT EXISTS "${failure_run}/${required_file}")
        message(FATAL_ERROR "failing contract omitted ${required_file}")
    endif()
endforeach()
file(READ "${failure_run}/audited.status" audited_failure_status)
string(STRIP "${audited_failure_status}" audited_failure_status)
if(NOT audited_failure_status STREQUAL "23")
    message(FATAL_ERROR "failing audited status was ${audited_failure_status}, expected 23")
endif()
file(READ "${failure_run}/contract.status" failure_contract_status)
string(STRIP "${failure_contract_status}" failure_contract_status)
if(NOT failure_contract_status STREQUAL "FAILED")
    message(FATAL_ERROR "failing child has contract status ${failure_contract_status}")
endif()
if(EXISTS "${failure_run}/contract.success")
    message(FATAL_ERROR "failing child left a success marker")
endif()

message(STATUS "Spin-It 2014 dependency isolation behavior passed: ${FIXTURE_ROOT}")
