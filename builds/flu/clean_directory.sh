#!/bin/bash

echo "Removing temporary directories."
rm -rf temp*
echo "Clearing out auspice/, processed/, prepared/"
rm auspice/*
rm processed/*
rm prepared/*
echo "Done cleaning directories, ready for clean builds."
