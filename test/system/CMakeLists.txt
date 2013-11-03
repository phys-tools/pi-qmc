set(SYSTEM_TEST_SCRIPTS
    README.md
    run_tests.py
    hatom/__init__.py
    hatom/pimc.xml
    hatom/test_hatom.py
    madelung/__init__.py
    madelung/paths.in
    madelung/pimc.xml
    madelung/test_nacl.py
    sho/__init__.py
    sho/pimc.xml
    sho/test_sho.py
)

foreach(_script ${SYSTEM_TEST_SCRIPTS})
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/${_script}
        ${CMAKE_CURRENT_BINARY_DIR}/${_script}
        COPYONLY
    )
endforeach()