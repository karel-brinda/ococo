OCOCO - Online Consensus Caller
===============================

Compilation
-----------

.. code-block:: bash

	git clone https://github.com/karel-brinda/ococo
	cd ococo
	cmake .
	make

Compilation for debugging mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Booost detailed logging and asserts.

.. code-block:: bash

	cmake -DDEBUGGING_MODE=OFF .
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

	Ococo: On-line consensus caller.:
	  -i [ --input ] arg      Input SAM/BAM file (- for standard input).
	  -f [ --fasta-ref ] arg  Initial FASTA reference (if not provided, sequence of
	                          N's is considered as the reference).
	  -s [ --stats ] arg      File with up-to-date statistics.
	  -S [ --strategy ] arg   Strategy for updates: majority / randomized. 
	                          [majority]
	  -m [ --mode ] arg       Mode: real-time / batch. [batch]
	  -q [ --min-MQ ] arg     Minimal mapping quality to increment a counter. [1]
	  -Q [ --min-BQ ] arg     Minimal base quality to increment a counter. [0]
	  -a [ --allow-amb ]      Ambiguous nucleotides can be called.
	  -w [ --ref-weight ] arg Initial counter value for nucleotides from reference.
	                          [2]
	  -v [ --vcf-cons ] arg   VCF file with updates of consensus.
	  -c [ --fasta-cons ] arg FASTA file with consensus, which is continuously 
	                          updated (WARNING: will be rewritten).
