/*
========================================================================================
    QUALITY CONTROL MODULES
========================================================================================
*/

process FASTQC {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
        'quay.io/biocontainers/fastqc:0.11.9--0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    if (meta.single_end) {
        """
        fastqc $args --threads $task.cpus $reads
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastqc: \$(fastqc --version | sed 's/FastQC v//')
        END_VERSIONS
        """
    } else {
        """
        fastqc $args --threads $task.cpus $reads
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastqc: \$(fastqc --version | sed 's/FastQC v//')
        END_VERSIONS
        """
    }
}

process MULTIQC {
    label 'process_single'

    conda (params.enable_conda ? "bioconda::multiqc=1.13" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.13--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0' }"

    input:
    path  multiqc_files, stageAs: "?/*"
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
    def logo = multiqc_logo ? "--cl-config 'custom_logo: \"$multiqc_logo\"'" : ''
    """
    multiqc \\
        -f \\
        $args \\
        $config \\
        $extra_config \\
        $logo \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}

process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_single'

    // Requires `pyyaml` which does not have a dedicated container but is in the MultiQC container
    conda (params.enable_conda ? "bioconda::multiqc=1.13" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.13--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0' }"

    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env python3

    import yaml
    import platform
    from textwrap import dedent

    def _make_versions_html(versions):
        html = [
            dedent(
                '''\\
                <style>
                #nf-core-versions tbody:nth-child(even) {
                    background-color: #f2f2f2;
                }
                </style>
                <table class="table" style="width:100%" id="nf-core-versions">
                    <thead>
                        <tr>
                            <th> Process Name </th>
                            <th> Software </th>
                            <th> Version  </th>
                        </tr>
                    </thead>
                '''
            )
        ]
        for process, tmp_versions in sorted(versions.items()):
            html.append("<tbody>")
            for i, (tool, version) in enumerate(sorted(tmp_versions.items())):
                html.append(
                    dedent(
                        f'''\\
                        <tr>
                            <td><samp>{process if (i == 0) else ""}</samp></td>
                            <td><samp>{tool}</samp></td>
                            <td><samp>{version}</samp></td>
                        </tr>
                        '''
                    )
                )
            html.append("</tbody>")
        html.append("</table>")
        return "".join(html)

    with open("$versions") as f:
        versions = yaml.safe_load(f)

    # aggregate versions by the module name
    module_versions = {}
    for process_name, process_versions in versions.items():
        module_name = process_name.split(":")[-1]
        if module_name not in module_versions:
            module_versions[module_name] = {}
        module_versions[module_name].update(process_versions)

    # dump to YAML
    with open("software_versions.yml", "w") as f:
        yaml.dump(module_versions, f, default_flow_style=False)

    with open("software_versions_mqc.yml", "w") as f:
        mqc_yml_out = {
            "id": "software_versions",
            "section_name": "Software Versions",
            "section_href": "https://github.com/lilei1/immune_neoantigen_pipeline",
            "plot_type": "html",
            "description": "are collected at run time from the software output.",
            "data": _make_versions_html(module_versions),
        }
        yaml.dump(mqc_yml_out, f, default_flow_style=False)

    with open("versions.yml", "w") as f:
        yaml.dump({"${task.process}": {"python": platform.python_version(), "yaml": yaml.__version__}}, f, default_flow_style=False)
    """
}
