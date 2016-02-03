OCOCO - Online Consensus Caller
===============================

Compilation
-----------

.. code-block:: bash

	git clone https://github.com/karel-brinda/ococo
	cd ococo
	cmake .
	make


Compilation options
~~~~~~~~~~~~~~~~~~~

All options:

* ``DEBUG`` - Booost detailed logging and asserts.
* ``VERBOSE_VCF`` - Verbose mode for VCF (unchanged bases also reported).
* ``OCOCO32`` - More accurate statistics (32bits / per position instead of 16bits).
* ``INSTALL_DEBUG_SCRIPTS`` - Install also auxiliary debugging scripts by ``make install``.

Example:

.. code-block:: bash

	cmake -DDEBUG=ON .
	make



Installation
------------

.. code-block:: bash
	
	make install


Tests
-----

.. code-block:: bash

	ctest --output-on-failure


Parameters
----------

.. code-block::

	Command-line parameters:
	  -i [ --input ] arg              Input SAM/BAM file (- for standard input).
	  -f [ --fasta-ref ] arg          Initial FASTA reference (if not provided, 
	                                  sequence of N's is considered as the 
	                                  reference).
	  -F [ --fasta-cons ] arg         FASTA file with consensus, which is 
	                                  continuously updated.
	  -s [ --stats-in ] arg           Input statistics.
	  -S [ --stats-out ] arg          Outputs statistics.
	  -v [ --vcf-cons ] arg           VCF file with updates of consensus.
	  -m [ --mode ] arg               Mode: real-time / batch. [batch]
	  -t [ --strategy ] arg           Strategy for updates: majority / stochastic. 
	                                  [stochastic]
	  -a [ --allow-amb ]              Allow updates to ambiguous nucleotides.
	  -q [ --min-MQ ] arg             Skip alignments with mapping quality smaller 
	                                  than INT. [1]
	  -Q [ --min-BQ ] arg             Skip bases with base quality smaller than 
	                                  INT. [13]
	  -w [ --ref-weight ] arg         Initial counter value for nucleotides from 
	                                  the reference. [2]
	  -c [ --min-coverage ] arg       Minimum coverage required for update. [2]
	  -M [ --majority-threshold ] arg Majority threshold. [0.6]