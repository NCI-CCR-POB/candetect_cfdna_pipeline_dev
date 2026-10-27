gc5mb="${1}"
blklist="${2}"
keep="${3}"

bedtools intersect -a $gc5mb -b $blklist -v -wa |
    grep -v _ |
    awk '{ if ($4 >= 0.3) print $0 }' > $keep


# This script is not in current use but may need to be used if there is not an appropriate keep_5mb.bed. To add this back to the smk place this rule:
#rule make_gc_map_bind:
#    container: frag_container,
#    input:
#        gc5mb = config["gc5mb"],
#        blklist = config["blklist"],
#    log: logdir + "/make_gc_map_bind.log",
#    output: refdir + "/keep_5mb.bed",
#    params:
#        script =frag_scriptdir + "/make_gc_map_bind.sh",
#    shell:
#        """
#        {params.script} \
#        {input.gc5mb} \
#        {input.blklist} \
#        {output} &> {log}
#        """