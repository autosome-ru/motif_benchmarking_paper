#!/usr/bin/env bash
MATRIX=$( echo $1 | grep -oPe '^MA\d+\.\d+' )
curl -qs "http://jaspar.genereg.net/api/v1/matrix/${MATRIX}/?format=json" | jq -r '.family|join(";")'
