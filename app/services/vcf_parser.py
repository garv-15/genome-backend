import os
import tempfile
from typing import List, Dict, Any
from cyvcf2 import VCF
from fastapi import UploadFile

# 0: HOM_REF (0/0), 1: HET (0/1), 2: HOM_ALT (1/1), 3: UNKNOWN (./.)
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
            read_depth = int(variant.gt_depths[0])
        elif 'DP' in variant.INFO:
            read_depth = int(variant.INFO['DP'])

        genotype_quality = None
        if variant.gt_quals is not None and len(variant.gt_quals) > 0 and variant.gt_quals[0] >= 0:
            genotype_quality = float(variant.gt_quals[0])

        chrom = str(variant.CHROM)
        if not chrom.startswith("chr") and chrom != "MT":
            chrom = f"chr{chrom}"

        for alt in variant.ALT:
            record = {
                "variant_id": f"{chrom}_{variant.POS}_{variant.REF}_{alt}",
                "chromosome": chrom,
                "position": int(variant.POS),
                "reference_allele": str(variant.REF),
                "alternate_allele": str(alt),
                "quality_score": round(float(variant.QUAL), 2) if variant.QUAL is not None else None,
                "read_depth": read_depth,
                "zygosity": zygosity,
                "genotype_quality": genotype_quality,
                "filter_status": str(variant.FILTER) if variant.FILTER else "PASS"
            }
            parsed_variants.append(record)

    vcf.close()
    return parsed_variants


async def parse_uploaded_vcf(upload_file: UploadFile) -> List[Dict[str, Any]]:
    suffix = ".vcf.gz" if upload_file.filename.endswith(".gz") else ".vcf"
    
    with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
        content = await upload_file.read()
        tmp.write(content)
        tmp_path = tmp.name

    try:
        variants = parse_vcf_file(tmp_path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    return variants