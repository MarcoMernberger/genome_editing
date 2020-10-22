graph
[
	node
    [
        id 0
        label "PFIE_/tmp/pytest-of-mernberger/pytest-66/test_set_primer_range0/dna/genome.fasta_prep_fasta"
        graphics
        [
            fill "#FFFFFF"
        ]
    ]
	node
    [
        id 1
        label "../../../../../../tmp/pytest-of-mernberger/pytest-66/test_set_primer_range0/dna/genome.fasta"
        graphics
        [
            fill "#5050AF"
        ]
    ]
	node
    [
        id 2
        label "PFIE_FBgenomeGTF_gtf_file"
        graphics
        [
            fill "#FFFFFF"
        ]
    ]
	node
    [
        id 3
        label "../../../../../../tmp/pytest-of-mernberger/pytest-66/test_set_primer_range0/cdna/cdna.fasta"
        graphics
        [
            fill "#5050AF"
        ]
    ]
	node
    [
        id 4
        label "../../../../../../tmp/pytest-of-mernberger/pytest-66/test_set_primer_range0/lookup/df_transcripts.msgpack"
        graphics
        [
            fill "#5050AF"
        ]
    ]
	node
    [
        id 5
        label "../../../../../../tmp/pytest-of-mernberger/pytest-66/test_set_primer_range0/lookup/df_genes.msgpack"
        graphics
        [
            fill "#5050AF"
        ]
    ]
	node
    [
        id 6
        label "../../../../../../tmp/pytest-of-mernberger/pytest-66/test_set_primer_range0/protein/pep.fasta"
        graphics
        [
            fill "#5050AF"
        ]
    ]
	node
    [
        id 7
        label "../../../../../../tmp/pytest-of-mernberger/pytest-66/test_set_primer_range0/lookup/df_proteins.msgpack"
        graphics
        [
            fill "#5050AF"
        ]
    ]
	node
    [
        id 8
        label "PysamCheck__ranges_by_amplicon"
        graphics
        [
            fill "#10FF10"
        ]
    ]
	node
    [
        id 9
        label "PysamCheck__lookup"
        graphics
        [
            fill "#10FF10"
        ]
    ]
	edge
    [
        source 0
        target 1
    ]
	edge
    [
        source 1
        target 3
    ]
	edge
    [
        source 4
        target 3
    ]
	edge
    [
        source 5
        target 4
    ]
	edge
    [
        source 2
        target 4
    ]
	edge
    [
        source 2
        target 5
    ]
	edge
    [
        source 7
        target 6
    ]
	edge
    [
        source 1
        target 6
    ]
	edge
    [
        source 2
        target 7
    ]
	edge
    [
        source 1
        target 8
    ]
	edge
    [
        source 9
        target 8
    ]
	edge
    [
        source 6
        target 8
    ]
	edge
    [
        source 3
        target 8
    ]
	edge
    [
        source 1
        target 9
    ]
	edge
    [
        source 6
        target 9
    ]
	edge
    [
        source 3
        target 9
    ]
]
