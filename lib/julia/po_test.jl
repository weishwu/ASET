using DataFrames
using CSV
using CodecZlib
using EstimatingEquationsRegression
using GLM
using LinearAlgebra
using Printf
using Statistics

function po_test(input_file::String) 
    println("Processing file: ", input_file)
    df = CSV.File(open(input_file); missingstring=["NA"]) |> DataFrame

    function prep(rout, eout, df)

        write(eout, @sprintf("%d distinct genes in initial file.\n", length(unique(df[:, :genes_exonic]))))

        f = r -> r.refAllele == r.PatAllele ? r.refCount : r.altCount
        df[:, :patCount] = [f(row) for row in eachrow(df)]

        f = r -> r.refAllele == r.MatAllele ? r.refCount : r.altCount
        df[:, :matCount] = [f(row) for row in eachrow(df)]

        # Phased data
        dp = filter(r -> !ismissing(r.PatAllele), df)
        dp = filter(r -> !ismissing(r.MatAllele), df)
        write(eout, @sprintf("%d distinct genes after dropping unphased data.\n", length(unique(dp[:, :genes_exonic]))))

        return dp
    end

    function build_pcs(fml, va, prefix, ratio, mx; u0 = nothing)
        fmls = apply_schema(fml, schema(fml, va))
        _, pred = modelcols(fmls, va)
        for j = 1:size(pred, 2)
            pred[:, j] .-= mean(pred[:, j])
        end

        if !isnothing(u0)
            # Project away from u0 if provided.
            u0, _, _ = svd(u0)
            pred = pred - u0 * u0' * pred
        end

        # Keep only the dominant factors
        u, s, _ = svd(pred)
        u = u[:, s.>=s[1]*ratio]
        if size(u, 2) > mx
            u = u[:, 1:mx]
        end

        # Standardize and create names
        pcn = String[]
        for k = 1:size(u, 2)
            u[:, k] .-= mean(u[:, k])
            u[:, k] ./= std(u[:, k])
            na = @sprintf("%s%d", prefix, k)
            push!(pcn, na)
            va[:, Symbol(na)] = u[:, k]
        end
        return va, u, pcn
    end

    function fitmodel(ky, v, rout, eout, block, po, pose, rap)

        # Create a dataframe of maternal read counts
        vm = v[:, [:RNAid, :variantID, :matCount, :refAllele, :MatAllele, :PatAllele]]
        vm = rename(vm, :matCount => :count)
        vm[:, :po] .= -1 / 2
        vm[:, :ra] = [r.refAllele == r.MatAllele ? 0.5 : -0.5 for r in eachrow(vm)]

        # Create a dataframe of paternal read counts
        vp = v[:, [:RNAid, :variantID, :patCount, :refAllele, :MatAllele, :PatAllele]]
        vp = rename(vp, :patCount => :count)
        vp[:, :po] .= 1 / 2
        vp[:, :ra] = [r.refAllele == r.PatAllele ? 0.5 : -0.5 for r in eachrow(vp)]

        # Merge the datasets
        va = vcat(vm, vp)
        va = sort(va, :RNAid)

        # Don't analyze genes with too few observations
        if size(va, 1) < 20
            d = size(va, 1)
            write(eout, @sprintf("[1] Skipping %s which has only %d observations.\n", string(ky), d))
            return block, po, pose, rap
        end

        # Base model without SNP effects
        f0 = @formula(count ~ po + ra)
        m0 = try
            gee(f0, va, va[:, :RNAid], Poisson(), IndependenceCor(), LogLink(); bccor = false)
        catch e
            write(
                eout,
                @sprintf("[2] Skipping %s because the base model did not fit.\n", string(ky))
            )
            return block, po, pose, rap
        end

        # If there is only one SNP, use the base model.
        if length(unique(va[:, :variantID])) == 1
            b0 = coef(m0)
            bse0 = sqrt.(diag(vcov(m0)))
            push!(block, ky)
            push!(po, b0[2])
            push!(pose, b0[2] / bse0[2])
            push!(rap, missing)
            write(eout, @sprintf("[3] %s has only one SNP so using base model.\n", string(ky)))
            return block, po, pose, rap
        end

        # Add SNP main effects
        va, mep, pcme = build_pcs(@formula(count ~ variantID), va, "PCME", 0.01, 10)
        tm = vcat([:po, :ra], Symbol.(pcme))
        f1 = Term(:count) ~ sum(Term.(tm))
        m1 = try
            gee(f1, va, va[:, :RNAid], Poisson(), IndependenceCor(), LogLink(); bccor = false)
        catch e
            write(eout, @sprintf("[4] Skipping %s because the main effects model did not fit.\n", string(ky)))
            return block, po, pose, rap
        end

        # Add SNP x ref/alt interactions, which capture genetic ASE
        va, ixp, pcix = build_pcs(@formula(count ~ variantID & ra), va, "PCIX", 0.01, 10; u0 = mep)
        tm = vcat(tm, Symbol.(pcix))
        f2 = Term(:count) ~ sum(Term.(tm))
        m2x = try
            gee(f2, va, va[:, :RNAid], Poisson(), IndependenceCor(), LogLink(), dofit=false; bccor=false)
        catch e
            write(eout, @sprintf("[5] Skipping %s because the interaction model did not fit.\n", string(ky)))
            return block, po, pose, rap
        end

        # Test SNP x refalt interactions versus base model with only SNP main effects
        st = scoretest(m2x.model, m1.model)

        m2 = try
            gee(f2, va, va[:, :RNAid], Poisson(), IndependenceCor(), LogLink(); bccor = false)
        catch e
            println(
                @sprintf(
                    "[6] Skipping %s because the interaction model did not fit.\n",
                    string(ky)
                )
            )
            return block, po, pose, rap
        end

        write(rout, @sprintf("Block=%s\n", string(ky)))
        write(
            rout,
            @sprintf("n=%d allele read counts (2 x number of SNPs x subjects)\n", nobs(m1))
        )
        write(rout, @sprintf("Score test for genetic ASE: p=%f\n", pvalue(st)))
        for (jj, mm) in enumerate([m0, m1, m2])
            ct = coeftable(mm)
            ct = split(string(ct), "\n")
            write(
                rout,
                [
                    "\nBase model\n",
                    "\nBase model + SNP main effects\n",
                    "\nBase model + SNP main effects + SNPxRefalt interactions\n",
                ][jj],
            )
            write(rout, join(ct, "\n"))
            write(rout, "\n")
        end
        write(rout, "\n")

        b2 = coef(m2)
        v = diag(vcov(m2))
        bse2 = sqrt.(v)

        push!(block, ky)
        push!(rap, pvalue(st))
        push!(po, b2[2])
        push!(pose, bse2[2])

        # Effect size for genetic ASE
        #b2 = b2[3:end]
        #bse2 = bse2[3:end]
        #b2 = b2[abs.(b2)./bse2.>2]

        return block, po, pose, rap
    end

    # Loop through the genes
    function main(dp, analysis_level, rout, eout)
        po, pose, rap = Union{Float64,Missing}[], Union{Float64,Missing}[], Union{Float64,Missing}[]
        if analysis_level == :gene
            block = Any[]
            for (k, v) in pairs(groupby(dp, :genes_exonic))
                block, po, pose, rap = fitmodel(k, v, rout, eout, block, po, pose, rap)
            end
            gene = first.(block)
            rslt = DataFrame(:genes_exonic => gene, :po => po, :po_se => pose, :Genetic_ASE_p => rap)
        elseif analysis_level == :exon
            block = Any[]
            for (k, v) in pairs(groupby(dp, [:genes_exonic, :exons_merged]))
                block, po, pose, rap = fitmodel(k, v, rout, eout, block, po, pose, rap)
            end
            gene = first.(block)
            exon = last.(block)
            rslt = DataFrame(:genes_exonic => gene, :exons_merged => exon, :po => po, :po_se => pose, :Genetic_ASE_p => rap)
        end

        rslt[:, :po_z] = rslt[:, :po] ./ rslt[:, :po_se]
        return rslt
    end


    for analysis_level in [:exon, :gene]

        rout = open("megpeg_models_$(analysis_level).txt", "w")
        eout = open("megpeg_models_$(analysis_level)_log.txt", "w")

        dp = prep(rout, eout, df)
        rslt = main(dp, analysis_level, rout, eout)
        CSV.write("megpeg_$(analysis_level).csv", rslt)

        close(rout)
        close(eout)
    end
end


po_test(ARGS[1])
