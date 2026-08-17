import os
import tempfile
import shutil
import numpy as np
from typing import List, Dict, Any
from cyvcf2 import VCF
from fastapi import UploadFile

ZYGOSITY_MAP = {
    0: "homozygous_ref",
    1: "heterozygous",
    2: "homozygous_alt",
    3: "unknown"
}

def parse_vcf_file(file_path: str) -> List[Dict[str, Any]]:
    parsed_variants = []
    vcf = VCF(file_path)

    for variant in vcf:
        gt_type = variant.gt_types[0] if len(variant.gt_types) > 0 else 3
        if gt_type == 0:
            continue

        zygosity = ZYGOSITY_MAP.get(gt_type, "unknown")

        read_depth = None
        if variant.gt_depths is not None and len(variant.gt_depths) > 0 and variant.gt_depths[0] >= 0:
            read_depth = float(variant.gt_depths[0])
        elif 'DP' in variant.INFO:
            read_depth = float(variant.INFO['DP'])

        genotype_quality = None
        if variant.gt_quals is not None and len(variant.gt_quals) > 0 and variant.gt_quals[0] >= 0:
            genotype_quality = float(variant.gt_quals[0])

        chrom = str(variant.CHROM)
        if not chrom.startswith("chr") and chrom not in ("MT", "M"):
            chrom = f"chr{chrom}"

        info = variant.INFO
        
        af = info.get('AF', np.nan)
        if isinstance(af, (tuple, list)):
            af = float(af[0]) if len(af) > 0 else np.nan
        elif af is not None:
            af = float(af)

        cadd_phred = float(info.get('CADD_PHRED', np.nan)) if 'CADD_PHRED' in info else np.nan
        revel_score = float(info.get('REVEL', np.nan)) if 'REVEL' in info else np.nan

        ref = str(variant.REF)
        var_type = variant.var_type

        for alt in variant.ALT:
            alt_str = str(alt)
            ref_len = len(ref)
            alt_len = len(alt_str)
            
            record = {
                "variant_id": f"{chrom}_{variant.POS}_{ref}_{alt_str}",
                "chromosome": chrom,
                "position": int(variant.POS),
                "reference_allele": ref,
                "alternate_allele": alt_str,
                "zygosity": zygosity,
                "filter_status": str(variant.FILTER) if variant.FILTER else "PASS",
                
                # ML Model Input Features
                "ref_len": ref_len,
                "alt_len": alt_len,
                "len_diff": alt_len - ref_len,
                "is_snp": 1 if var_type == 'snp' else 0,
                "is_indel": 1 if var_type == 'indel' else 0,
                "qual": float(variant.QUAL) if variant.QUAL is not None else np.nan,
                "is_pass": 1 if (variant.FILTER is None or variant.FILTER == 'PASS') else 0,
                "allele_freq": af,
                "depth": read_depth if read_depth is not None else np.nan,
                "cadd_phred": cadd_phred,
                "revel_score": revel_score
            }
            parsed_variants.append(record)

    vcf.close()
    return parsed_variants


async def parse_uploaded_vcf(upload_file: UploadFile) -> List[Dict[str, Any]]:
    suffix = ".vcf.gz" if upload_file.filename.endswith(".gz") else ".vcf"
    
    with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
        shutil.copyfileobj(upload_file.file, tmp)
        tmp_path = tmp.name

    try:
        variants = parse_vcf_file(tmp_path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    return variants