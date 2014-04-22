  # A series of examples to demonstrate SWAMP's functions.
  # To run example1 type 'make example1'
  # To run example2 type 'make example2' etc.
  # Run 'make clean' to remove the generated masked files.
  # Further commentary on the examples can be found in README.md

example1:
  # Basic example of running SWAMP on a single alignment, with a large
  # window and low substitution threshold. These paramters are
  # purposefully over-aggressive to show the filtered alignment.

	python SWAMP.py --print-alignment example_dataset/data/44/44.phy

	python SWAMP.py -i example_dataset/data/44/ \
	-b example_dataset/branchcodes.txt -t 1 -w 30

	python SWAMP.py --print-alignment example_dataset/data/44/44_masked.phy


example2:
  # Run SWAMP with more reasonable parameters, summarizing the masked
  # alignment.

	python SWAMP.py -i example_dataset/data/44/ \
	-b example_dataset/branchcodes.txt -t 2 -w 20

	python SWAMP.py --print-alignment example_dataset/data/44/44_masked.phy

example3:
  # Run with interscan=true, which recursively merges filtered regions
  # to increase the filtering stringency.

	python SWAMP.py -i example_dataset/data/44/ \
	-b example_dataset/branchcodes.txt -t 2 -w 20 --interscan

	python SWAMP.py --print-alignment example_dataset/data/44/44_masked.phy


example4:
  # Run with very stringent filtering on just the human sequence.
	python SWAMP.py -i example_dataset/data/44/ \
	-b example_dataset/branchcodes_onlyhuman.txt -t 2 -w 50 --interscan

  # No codons were masked.
	python SWAMP.py --print-alignment example_dataset/data/44/44_masked.phy

example5:
  # Run on the entire directory within example_dataset/data/
  # The branchcodes_nohomo.txt file in this example will not mask the 
  # homo sequences as they are not listed.
	python SWAMP.py -i example_dataset/data/ \
	-b example_dataset/branchcodes.txt -t 5 -w 15 --interscan

clean:
  # Clean out all masked files from example dataset directories.
	rm -v example_dataset/data/*/*_masked.phy

# Just tells make that these tasks should be run even though they have
# no source files.
.PHONY: example1 example2 example3 example4 example5 clean
