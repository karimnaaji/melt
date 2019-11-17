file(READ "${in}" file_content)
string(REGEX REPLACE "([\n]*[ \t]*(#)[^\n]*)" "" content "${file_content}")
string(REGEX REPLACE "[\n]+" "\n" content "${content}")
file(WRITE "${out}" "#pragma once\n\nconst static char s_${name}[] = R\"(\n${content})\";\n")
