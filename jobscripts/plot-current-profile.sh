#!/bin/bash
# thin wrapper: the current profiles are set up and run by gimic-profile (installed with GIMIC)
exec gimic-profile --plot "${1:-.}" "${@:2}"
