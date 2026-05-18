# fastder input preparation, build and run rules.
#
# Included by the main Snakefile after the path constants and the parameter
# grid block are defined.


# 7b. Extract and organise all inputs that fastder needs into a flat directory.
# The set of inputs depends on the backend: with monorail we read from the
# unify/pump output dirs; with monorail_light we read from the LIGHT_DIR
# scratch produced by the ml_* rules above.
def _extract_inputs(wc):
    if BACKEND == "recount3":
        # wc.scenario is a recount3 sample group. Inputs are the group's
        # lean MM/RR (built by recount3_group_junctions) and the downloaded
        # per-sample BigWigs for that group.
        if STRANDED:
            raise ValueError(
                "fastder.stranded is not supported with the recount3 backend: "
                "recount3 hosts unstranded coverage only."
            )
        group_samples = R3_GROUPS[wc.scenario]
        return {
            "rr": op.join(R3_DIR, wc.scenario, "junctions.ALL.RR"),
            "mm": op.join(R3_DIR, wc.scenario, "junctions.ALL.MM"),
            "samples_tsv": op.join(R3_DIR, wc.scenario, "junctions.ALL.samples.tsv"),
            "bws": [op.join(R3_DIR, "bw", f"{s}.all.bw") for s in group_samples],
            # The reference annotation is the gffcompare truth set; depend on
            # the download so it is present before extract_fastder_inputs runs.
            "reference_gtf": REF_GTF,
        }
    if BACKEND == "monorail_light":
        light_scn = op.join(LIGHT_DIR, wc.scenario)
        result = {
            "rr": op.join(light_scn, "junctions.ALL.RR"),
            "mm": op.join(light_scn, "junctions.ALL.MM"),
            "samples_tsv": op.join(light_scn, "samples.tsv"),
        }
        if STRANDED:
            result["bws"] = expand(
                [op.join(light_scn, "{sample}.plus.bw"),
                 op.join(light_scn, "{sample}.minus.bw")],
                sample=PUMP_SAMPLES,
            )
        else:
            result["bws"] = expand(op.join(light_scn, "{sample}.all.bw"), sample=PUMP_SAMPLES)
        # The per-sample truth GFFs exist only for ASimulatoR input. Local
        # FASTQ input on this backend has no simulated truth.
        if PUMP_SOURCE == "asimulator":
            result["sample_gffs"] = expand(
                op.join(DATA_DIR, "asim", "{sample}", "{scenario}", "splicing_variants.gff3"),
                sample=PUMP_SAMPLES, scenario=[wc.scenario],
            )
        return result
    # monorail (heavy) backend (scenario unused; legacy path)
    result = {
        "unify_dir": op.join(DATA_DIR, "unify"),
        "pump_dirs": expand(op.join(DATA_DIR, "pump", "{sample}"), sample=PUMP_SAMPLES),
    }
    if STRANDED:
        result["stranded_bws"] = expand(
            [op.join(DATA_DIR, "stranded_bigwigs", "{sample}.plus.bw"),
             op.join(DATA_DIR, "stranded_bigwigs", "{sample}.minus.bw")],
            sample=PUMP_SAMPLES,
        )
    # SRA and local input are graded against the reference annotation. Depend
    # on the download so the annotation is present before the run body needs
    # it (ASimulatoR input gets it through run_asimulator instead).
    if PUMP_SOURCE in ("sra", "local"):
        result["reference_gtf"] = REF_GTF
    return result


rule extract_fastder_inputs:
    input:
        unpack(_extract_inputs),
    output:
        touch(op.join(FASTDER_DIR, "{scenario}", "extract.DONE"))
    benchmark:
        op.join(BENCH_DIR, "extract_fastder_inputs_{scenario}.tsv")
    log:
        op.join(LOG_DIR, "extract_fastder_inputs_{scenario}.log")
    params:
        fastder_dir=lambda wc: op.join(FASTDER_DIR, wc.scenario),
        pump_dir=op.join(DATA_DIR, "pump"),
        light_dir=lambda wc: op.join(LIGHT_DIR, wc.scenario),
        scenario=lambda wc: wc.scenario,
        scenario_samples=lambda wc: SAMPLES_BY_SCENARIO[wc.scenario],
        asim_dir=op.join(DATA_DIR, "asim"),
        project_name=config["monorail"]["project_name"],
        ref_version=config["monorail"]["ref_version"],
        stranded=STRANDED,
        stranded_bw_dir=op.join(DATA_DIR, "stranded_bigwigs"),
        samples_tsv_script=op.join(WORKFLOW_DIR, "scripts", "create_bigwig_list.sh"),
        backend=BACKEND,
    run:
        import shutil, gzip

        fdir = Path(params.fastder_dir)
        fdir.mkdir(parents=True, exist_ok=True)

        if params.backend == "recount3":
            # 1. Copy the group's lean MM/RR straight in (already subset to
            # the group's samples and filtered to the analysed chromosomes).
            shutil.copy2(input.rr, fdir / "junctions.ALL.RR")
            shutil.copy2(input.mm, fdir / "junctions.ALL.MM")

            # 2. Copy the group's downloaded BigWigs. recount3_fetch_bigwig
            # already names them {sample}.all.bw, which is what fastder wants.
            for bw in input.bws:
                shutil.copy2(bw, fdir / Path(bw).name)

            # 3. Generate the BigWig-list CSV from the group samples.tsv.
            bigwig_csv = fdir / "recount3_BigWig_list.csv"
            shell(
                "bash {params.samples_tsv_script}"
                f" {input.samples_tsv}"
                f" {bigwig_csv}"
                " unstranded"
                " >> {log} 2>&1"
            )
        elif params.backend == "monorail_light":
            # 1. Copy lean MM/RR straight in (already gunzipped, already filtered)
            shutil.copy2(input.rr, fdir / "junctions.ALL.RR")
            shutil.copy2(input.mm, fdir / "junctions.ALL.MM")

            # 2. Copy BigWigs (the ml_bam_to_bigwig rule already named them
            # {sample}.{strand}.bw / {sample}.all.bw to match what fastder expects)
            ldir = Path(params.light_dir)
            for sample in PUMP_SAMPLES:
                if params.stranded:
                    for strand in ("plus", "minus"):
                        shutil.copy2(ldir / f"{sample}.{strand}.bw",
                                     fdir / f"{sample}.{strand}.bw")
                else:
                    shutil.copy2(ldir / f"{sample}.all.bw", fdir / f"{sample}.all.bw")

            # 3. Generate BigWig-list CSV from the lean samples.tsv
            bigwig_csv = fdir / "recount3_BigWig_list.csv"
            bw_mode = "stranded" if params.stranded else "unstranded"
            shell(
                "bash {params.samples_tsv_script}"
                f" {input.samples_tsv}"
                f" {bigwig_csv}"
                f" {bw_mode}"
                " >> {log} 2>&1"
            )
        else:
            # 1. Gunzip MM and RR files from unify output
            project = params.project_name
            jct_base = Path(input.unify_dir) / "junction_counts_per_study"
            mm_gz = rr_gz = None
            for gz in jct_base.rglob(f"*.ALL.MM.gz"):
                mm_gz = gz
            for gz in jct_base.rglob(f"*.ALL.RR.gz"):
                rr_gz = gz

            if not mm_gz or not rr_gz:
                raise FileNotFoundError(
                    f"Could not find ALL.MM.gz / ALL.RR.gz under {jct_base}"
                )

            for gz_path, suffix in [(mm_gz, ".ALL.MM"), (rr_gz, ".ALL.RR")]:
                out_name = gz_path.name.removesuffix(".gz")
                with gzip.open(gz_path, "rb") as f_in, open(fdir / out_name, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)

            # 2. Copy BigWig files into fastder directory
            if params.stranded:
                bw_dir = Path(params.stranded_bw_dir)
                for sample in PUMP_SAMPLES:
                    for strand in ("plus", "minus"):
                        bw_src = bw_dir / f"{sample}.{strand}.bw"
                        if not bw_src.exists():
                            raise FileNotFoundError(f"Stranded BigWig not found: {bw_src}")
                        shutil.copy2(bw_src, fdir / f"{sample}.{strand}.bw")
            else:
                for sample in PUMP_SAMPLES:
                    att_dir = Path(params.pump_dir) / sample / "output" / f"{sample}_att0"
                    if PUMP_SOURCE in ("asimulator", "local"):
                        # Both pump local files and label the study "local".
                        study = project
                    else:
                        study = config["monorail"]["sra_samples"][sample]["study_acc"]
                    bw_name = f"{sample}!{study}!{params.ref_version}!local.all.bw"
                    bw_src = att_dir / bw_name
                    if not bw_src.exists():
                        raise FileNotFoundError(f"BigWig not found: {bw_src}")
                    shutil.copy2(bw_src, fdir / f"{sample}.all.bw")

            # 3. Generate BigWig-list CSV from unify samples.tsv
            samples_tsv = Path(input.unify_dir) / "samples.tsv"
            bigwig_csv = fdir / "recount3_BigWig_list.csv"
            bw_mode = "stranded" if params.stranded else "unstranded"
            shell(
                "bash {params.samples_tsv_script}"
                f" {samples_tsv}"
                f" {bigwig_csv}"
                f" {bw_mode}"
                " >> {log} 2>&1"
            )

        # 4. Copy label / reference GFF for gffcompare, one per sample that
        # this scenario is evaluated against. With the ASimulatoR backends the
        # truth set is the per-scenario splicing_variants.gff3, so the
        # variant_only truth contains only the alternative isoforms. With the
        # recount3 backend the truth is the reference annotation and there is
        # one pseudo-sample, "reference".
        if HAS_SIM_TRUTH:
            for sample in params.scenario_samples:
                gff_src = (Path(str(params.asim_dir)) / sample
                           / params.scenario / "splicing_variants.gff3")
                if not gff_src.exists():
                    gff_src = Path(str(params.asim_dir)) / sample / "splicing_variants.gff3"
                if gff_src.exists():
                    shutil.copy2(gff_src, fdir / f"{sample}_label.gff3")
        else:
            ref_path = Path(str(REF_ANNOTATION))
            if not ref_path.is_absolute():
                ref_path = Path(str(WORKFLOW_DIR)) / ref_path
            if not ref_path.exists():
                raise FileNotFoundError(f"reference annotation not found: {ref_path}")
            for sample in params.scenario_samples:
                shutil.copy2(ref_path, fdir / f"{sample}_label{ref_path.suffix}")


# 8. Match chr prefix convention for files used by gffcompare
rule match_chr_prefix:
    input:
        op.join(FASTDER_DIR, "{scenario}", "extract.DONE")
    output:
        touch(op.join(FASTDER_DIR, "{scenario}", "match_chr_prefix.DONE"))
    benchmark:
        op.join(BENCH_DIR, "match_chr_prefix_{scenario}.tsv")
    log:
        op.join(LOG_DIR, "match_chr_prefix_{scenario}.log")
    params:
        fastder_dir=lambda wc: op.join(FASTDER_DIR, wc.scenario),
    conda:
        "../envs/base.yaml"
    shell:
        """
        for file in {params.fastder_dir}/*_label.gff3 {params.fastder_dir}/*_label.gtf; do
            [ -f "$file" ] || continue
            tmp="${{file}}.tmp"
            awk 'BEGIN{{FS=OFS="\\t"}}
                /^#/ {{ print; next }}
                $1 !~ /^chr/ {{ $1 = "chr" $1 }}
                {{ print }}
            ' "$file" > "$tmp"
            mv "$tmp" "$file"
        done > {log} 2>&1
        """


# Legacy rule, not wired into the active DAG.
#
# This rule used to convert each *.bw file in data/fastder/ into a *.bedGraph
# so that older fastder builds without libBigWig support could read the
# coverage. The current build_fastder configures fastder with
# -DFASTDER_USE_LIBBIGWIG=ON and run_fastder symlinks the .bw files directly,
# so this conversion is unnecessary on the active code path. The rule is kept
# in place for two reasons:
#   1. If libBigWig support is intentionally disabled (e.g. for a packaging
#      build that does not pull libBigWig), wiring this rule back into
#      run_fastder's input list and switching the symlink glob to *.bedGraph
#      restores the old BedGraph-fed pipeline with no further changes.
#   2. It documents the historical conversion step that produced the
#      BedGraph fixtures used in earlier benchmarks.
# To re-enable: add
# `bedgraph_done=op.join(FASTDER_DIR, "{scenario}", "convert_bedgraph.DONE")`
# to run_fastder's input clause and change its symlink loop to *.bedGraph.
# Paths carry the {scenario} subdirectory so the rule composes with the
# scenario-aware layout introduced by the make_scenario rule.
rule bigwig_to_bedgraph:
    input:
        op.join(FASTDER_DIR, "{scenario}", "extract.DONE")
    output:
        touch(op.join(FASTDER_DIR, "{scenario}", "convert_bedgraph.DONE"))
    benchmark:
        op.join(BENCH_DIR, "bigwig_to_bedgraph_{scenario}.tsv")
    log:
        op.join(LOG_DIR, "bigwig_to_bedgraph_{scenario}.log")
    params:
        fastder_dir=lambda wc: op.join(FASTDER_DIR, wc.scenario),
    conda:
        "../envs/ucsc_tools.yaml"
    shell:
        """
        for bw in {params.fastder_dir}/*.bw; do
            [ -f "$bw" ] || continue
            out="${{bw%.bw}}.bedGraph"
            bigWigToBedGraph "$bw" "$out"
        done > {log} 2>&1
        """


# Build fastder from source. The build is configured with
# FASTDER_USE_LIBBIGWIG=ON so the binary reads .bw files directly without an
# intermediate BedGraph conversion step. libBigWig is fetched at configure
# time and needs zlib and libcurl (provided by envs/fastder_build.yaml).
rule build_fastder:
    output:
        op.join(FASTDER_BUILD_DIR, "fastder")
    benchmark:
        op.join(BENCH_DIR, "build_fastder.tsv")
    log:
        op.join(LOG_DIR, "build_fastder.log")
    params:
        fastder_src=op.join(WORKFLOW_DIR, "external", "fastder"),
        build_dir=str(FASTDER_BUILD_DIR),
        ncores=config["cores"],
    conda:
        "../envs/fastder_build.yaml"
    shell:
        """
        if [ -n "{params.build_dir}" ] && [ -d "{params.build_dir}" ]; then
            rm -rf "{params.build_dir}"
        fi
        mkdir -p {params.build_dir}
        cd {params.build_dir}
        cmake -DFASTDER_USE_LIBBIGWIG=ON {params.fastder_src} > {log} 2>&1
        cmake --build . --target fastder -j {params.ncores} >> {log} 2>&1
        """


# 10. Run fastder for each parameter combination.
# Each run gets its own working directory (data/fastder/runs/{param_id}/) with
# symlinks to the shared inputs. This allows parallel execution without races
# on the FASTDER_RESULT_*.gtf output filenames. Symlinks include .bw files
# directly; fastder reads them via libBigWig.
rule run_fastder:
    input:
        fastder_exe=op.join(FASTDER_BUILD_DIR, "fastder"),
        extract_done=op.join(FASTDER_DIR, "{scenario}", "extract.DONE"),
    output:
        gtf_path=op.join(FASTDER_DIR, "{scenario}", "run_fastder_{param_id}.gtf_path"),
        done=touch(op.join(FASTDER_DIR, "{scenario}", "run_fastder_{param_id}.DONE")),
    benchmark:
        op.join(BENCH_DIR, "run_fastder", "{scenario}_{param_id}.tsv")
    log:
        op.join(LOG_DIR, "run_fastder", "{scenario}_{param_id}.log")
    params:
        fastder_dir=lambda wc: op.join(FASTDER_DIR, wc.scenario),
        run_dir=lambda wc: op.join(FASTDER_DIR, wc.scenario, "runs", wc.param_id),
        ncores=config["cores"],
        fastder_args=lambda wc: PARAM_CLI_ARGS[wc.param_id],
        stranded_arg="--stranded" if STRANDED else "",
        chr_args=(
            "--chr " + " ".join(str(c) for c in FASTDER_CFG["chromosomes"])
            if FASTDER_CFG.get("chromosomes")
            else ""
        ),
    conda:
        "../envs/fastder_build.yaml"
    shell:
        """
        mkdir -p {params.run_dir}
        # Symlink all shared inputs into the private run directory.
        # Includes .bw (parsed by fastder via libBigWig), MM, RR, and the CSV.
        for f in {params.fastder_dir}/*.bw \
                  {params.fastder_dir}/*.MM \
                  {params.fastder_dir}/*.RR \
                  {params.fastder_dir}/*.csv; do
            [ -e "$f" ] || continue
            ln -sf "$f" {params.run_dir}/
        done
        {input.fastder_exe} \
            --dir {params.run_dir} \
            {params.stranded_arg} \
            {params.fastder_args} \
            --cores {params.ncores} \
            {params.chr_args} \
            > {log} 2>&1
        gtf=$(ls {params.run_dir}/FASTDER_RESULT_*.gtf 2>/dev/null)
        if [ -z "$gtf" ]; then
            echo "ERROR: fastder did not produce a FASTDER_RESULT_*.gtf" >&2; exit 1
        fi
        echo "$gtf" > {output.gtf_path}
        """


# 11b. Per-tool runners. Each rule produces a GTF at
# data/tools/{tool}/{scenario}/{param_id}/output.gtf so run_gffcompare and
# eval_fuzzy_metrics can compare the methods on the same simulated truth.

# fastder: re-export the run_fastder GTF at the standardised path.
rule link_fastder_gtf:
    input:
        gtf_path=op.join(FASTDER_DIR, "{scenario}", "run_fastder_{param_id}.gtf_path"),
    output:
        gtf=op.join(DATA_DIR, "tools", "fastder", "{scenario}", "{param_id}", "output.gtf"),
    benchmark:
        op.join(BENCH_DIR, "link_fastder_gtf", "{scenario}_{param_id}.tsv")
    log:
        op.join(LOG_DIR, "link_fastder_gtf", "{scenario}_{param_id}.log"),
    conda:
        "../envs/base.yaml"
    shell:
        """
        mkdir -p $(dirname {output.gtf})
        cp -f $(cat {input.gtf_path}) {output.gtf} 2> {log}
        """
