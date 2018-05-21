#!/bin/bash

echo "Removing temporary directories."
rm -rf temp*
echo "Clearing out prepared/, processed/, auspice/"
rm -rf prepared/*
rm -rf processed/*
rm auspice/*
echo "Done cleaning directories, ready for clean builds."
