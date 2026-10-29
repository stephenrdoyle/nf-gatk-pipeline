process SOFTWARE_VERSIONS {
    /*
     * Record the software this run used, for methods sections and for
     * checking whether two runs are comparable:
     * results/pipeline_info/software_versions.txt
     *
     * Runs under the same process.beforeScript as the GATK/samtools/bcftools
     * steps and loads nothing itself, so it reports the environment those
     * steps actually get rather than a second, hand-maintained list. Every
     * loaded module is listed, with its Singularity image for the
     * wrapper-based ones. MultiQC loads its own module; its version is in
     * the MultiQC report.
     *
     * cache false: the module environment is not part of Nextflow's cache
     * key, so a cached copy could describe an older environment. It takes
     * seconds, so it simply reruns every time.
     */
    tag "versions"
    label 'process_low'
    cache false

    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    val pipeline_info   // pipeline, git, Nextflow and command-line lines from main.nf

    output:
    path "software_versions.txt", emit: versions

    script:
    """
    {
        echo "# Software versions"
        echo "date: \$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        cat << 'PIPELINE_INFO'
${pipeline_info}
PIPELINE_INFO

        echo
        echo "## Tools"
        echo "gatk:     \$(gatk --version 2>&1 | grep -E '^(The Genome Analysis Toolkit|HTSJDK|Picard)' | paste -sd ';' -)"
        echo "samtools: \$(samtools --version 2>&1 | head -2 | paste -sd ';' -)"
        echo "bcftools: \$(bcftools --version 2>&1 | head -2 | paste -sd ';' -)"
        echo "python3:  \$(python3 --version 2>&1)"

        echo
        echo "## Modules (process.beforeScript)"
        # beforeScript defines 'module'; source it again only if this shell lacks it
        if ! type module >/dev/null 2>&1 && [ -f /etc/profile.d/modules.sh ]; then
            source /etc/profile.d/modules.sh
        fi
        for m in \$(module -t list 2>&1 | grep -v ':\$'); do
            img=\$(module show "\$m" 2>&1 | grep -o 'Singularity image location: [^}]*' | sed 's/.*: //')
            echo "\$m\${img:+    image: \$img}"
        done
    } > software_versions.txt
    """
}
