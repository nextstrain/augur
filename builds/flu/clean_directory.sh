#!/bin/bash

echo "Removing temporary directories."
rm -rf temp*
echo "Clearing out prepared/, processed/, auspice/"
rm prepared/*
rm processed/*
rm processed/recurring_mutations/*
rm processed/rising_mutations/*
rm auspice/*
echo "Done cleaning directories, ready for clean builds."
