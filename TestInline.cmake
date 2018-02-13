# Inspired from /usr/share/autoconf/autoconf/c.m4
FOREACH(KEYWORD "inline" "__inline__" "__inline")
    IF(NOT DEFINED C_INLINE)
        TRY_COMPILE(C_HAS_${KEYWORD} "${gmmem_BINARY_DIR}"
                "${gmmem_SOURCE_DIR}/test_inline.c"
                COMPILE_DEFINITIONS "-Dinline=${KEYWORD}")
        IF(C_HAS_${KEYWORD})
            SET(C_INLINE TRUE)
            ADD_DEFINITIONS("-Dinline=${KEYWORD}")
        ENDIF(C_HAS_${KEYWORD})
    ENDIF(NOT DEFINED C_INLINE)
ENDFOREACH(KEYWORD)
IF(NOT DEFINED C_INLINE)
    ADD_DEFINITIONS("-Dinline=")
ENDIF(NOT DEFINED C_INLINE)
