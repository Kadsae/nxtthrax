# snakefile for generating nextstrain output of anthracis phylogenetic analysis
# author: florian himmelbauer
# date created: 15-09-23
# date last modified: 28-01-25

rule all:
	input:
		"08_augur_out_export2auspice.json"


rule files:
	input:
		vcf = "01_parsnp_out/parsnp.vcf",
		reference = "data/reference/Ames-Ancestor.fna",
		annotation = "data/reference/Ames-Ancestor.gff3",
		lat_longs = "config/lat_longs.tsv",
		clades = "config/clades.tsv",
		metadata = "config/metadata.tsv",
		colors = "config/colors.tsv",


rule tree:
	message:
		"""
		Building tree from vcf
		"""
	input:
		alignment = rules.files.input.vcf,
		reference = rules.files.input.reference
	output:
		tree = "02_augur_out_tree_raw.nwk"
	shell:
		"""
		augur tree \
			--alignment {input.alignment} \
			--vcf-reference {input.reference} \
			--output {output.tree}
		"""


rule refine:
	message:
		"""
		Refine tree: Fix branch lenghts and get a time-resolved tree
		"""
	input:
		tree = rules.tree.output.tree,
		reference = rules.files.input.reference,
		alignment = rules.files.input.vcf,
		metadata = rules.files.input.metadata
	output:
		tree = "03_augur_out_tree_refined.nwk",
		node = "03_augur_out_tree_refined_branch_lengths.json"
	params:
		coalescent = "opt"
	shell:
		"""
		augur refine \
			--tree {input.tree} \
			--alignment {input.alignment} \
			--vcf-reference {input.reference} \
			--timetree \
			--root min_dev \
			--metadata {input.metadata} \
			--coalescent {params.coalescent} \
			--output-tree {output.tree} \
			--output-node-data {output.node}
		"""


rule ancestral:
	message:
		"""
		Reconstructing ancestral sequences
		"""
	input:
		tree = rules.refine.output.tree,
		alignment = rules.files.input.vcf,
		reference = rules.files.input.reference
	output:
		node = "04_augur_out_nt_muts.json",
		vcf = "04_augur_out_nt_muts.vcf"
	params:
		inference = "joint"
	shell:
		"""
		augur ancestral \
			--tree {input.tree} \
			--alignment {input.alignment} \
			--vcf-reference {input.reference} \
			--inference {params.inference} \
			--output-node-data {output.node} \
			--output-vcf {output.vcf} 
		"""


#rule translate:
#	message:
#		"""
#		Naming CDS regions to show nucleotide diversity of reference genome
#		"""
#	input:
#		tree = rules.refine.output.tree,
#		ancestral = rules.ancestral.output.vcf,
#		reference = rules.files.input.reference,
#		annotation = rules.files.input.annotation
#	output:
#		node = "05_augur_out_aa_muts.json"
#	shell:
#		"""
#		augur translate \
#			--tree {input.tree} \
#			--ancestral-sequences {input.ancestral} \
#			--vcf-reference {input.reference} \
#			--reference-sequence {input.annotation} \
#			--output-node-data {output.node}
#		"""


rule traits:
	message:
		"""
		Reconstruct ancestral traits
		"""
	input:
		tree = rules.refine.output.tree,
		metadata = rules.files.input.metadata
	output:
		node = "06_augur_out_traits.json"
	params:
		columns = "country location"
	shell:
		"""
		augur traits \
			--tree {input.tree} \
			--metadata {input.metadata} \
			--columns {params.columns} \
			--output-node-data {output.node}
		"""


rule clades:
	message:
		"""
		Identify specified clades
		"""
	input:
		tree = rules.refine.output.tree,
		clades = rules.files.input.clades,
		mutations = rules.ancestral.output.node
	output:
		node = "07_augur_out_clades.json"
	shell:
		"""
		augur clades \
			--tree {input.tree} \
			--clades {input.clades} \
			--mutations {input.mutations} \
			--output-node-data {output.node}
		"""


rule export:
	message:
		"""
		Export the results
		"""
	input:
		tree = rules.refine.output.tree,
		metadata = rules.files.input.metadata,
		node_tree = rules.refine.output.node,
		node_nt_muts = rules.ancestral.output.node,
		#node_aa_muts = rules.translate.output.node,
		node_traits = rules.traits.output.node,
		node_clades = rules.clades.output.node,
		lat_longs = rules.files.input.lat_longs,
		colors = rules.files.input.colors
	output:
		export = "08_augur_out_export2auspice.json"
	params:
		color = "host isolation_source AssemblyStatus checkm_completeness checkm_contamination data_origin sample_type"
	shell:
		"""
		augur export v2 \
			--tree {input.tree} \
			--metadata {input.metadata} \
			--color-by-metadata {params.color} \
			--colors {input.colors} \
			--node-data {input.node_tree} {input.node_nt_muts} {input.node_traits} {input.node_clades} \
			--geo-resolutions country location \
			--lat-longs {input.lat_longs} \
			--output {output.export}
		"""
