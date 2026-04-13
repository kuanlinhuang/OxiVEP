# Deploying fastVEP on Hetzner Cloud

Complete guide to hosting fastVEP as a public web service with full human genome annotation (GRCh38) and all supplementary databases.

## Cost Overview

| Component | Monthly Cost |
|-----------|-------------|
| Hetzner CCX23 (4 vCPU, 16GB RAM, 80GB SSD) | ~$30 |
| Cloudflare DNS + CDN (free tier) | $0 |
| **Total** | **~$30/month** |

Upgrade to CCX33 (8 vCPU, 32GB RAM, 160GB SSD, ~$55/month) if you want headroom for multiple genomes or heavy traffic.

## 1. Provision the Server

### Create account and server

1. Sign up at [hetzner.com/cloud](https://www.hetzner.com/cloud)
2. Create a new project
3. Add an SSH key (Settings → SSH Keys)
4. Create a server:
   - **Location**: Ashburn (us-east) or Falkenstein (eu-central) — pick closest to your users
   - **Image**: Ubuntu 24.04
   - **Type**: CCX23 (Dedicated, 4 vCPU AMD, 16GB RAM, 80GB SSD)
   - **SSH key**: select yours
   - **Name**: `fastvep`

### Initial setup

```bash
ssh root@YOUR_SERVER_IP

# System updates
apt update && apt upgrade -y

# Create service user
useradd -m -s /bin/bash fastvep
mkdir -p /opt/fastvep/data /opt/fastvep/bin
chown -R fastvep:fastvep /opt/fastvep

# Install nginx
apt install -y nginx certbot python3-certbot-nginx

# Firewall
ufw allow OpenSSH
ufw allow 'Nginx Full'
ufw enable
```

## 2. Build and Upload the Binary

On your development machine:

```bash
# Cross-compile for Linux (if building on Mac)
# Option A: Use cross (recommended)
cargo install cross
cross build --release --target x86_64-unknown-linux-gnu -p fastvep-web

# Option B: Build on the server directly
# (requires Rust on server: curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh)

# Upload binary
scp target/x86_64-unknown-linux-gnu/release/fastvep-web root@YOUR_SERVER_IP:/opt/fastvep/bin/
ssh root@YOUR_SERVER_IP 'chmod +x /opt/fastvep/bin/fastvep-web'
```

## 3. Download Genome Data

fastVEP supports multiple genomes via the `--data-dir` flag. Each subdirectory is an independent genome (different species or different human reference builds). The web UI shows them in a dropdown labeled with the build/source so users pick unambiguously.

### Directory layout

```
/opt/fastvep/data/
├── hg38_ensembl_115/              # ← each subdir = one genome
│   ├── Homo_sapiens.GRCh38.115.gff3
│   ├── Homo_sapiens.GRCh38.dna.primary_assembly.fa
│   ├── Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
│   └── sa/
│       ├── clinvar.osa
│       ├── gnomad.osa
│       └── ...
├── hg19_ensembl_113/
│   ├── Homo_sapiens.GRCh37.113.gff3
│   └── ...
├── mouse_grcm39_ensembl_115/
│   ├── Mus_musculus.GRCm39.115.gff3
│   └── ...
└── ...
```

The dropdown in the web UI will show each directory name (so **name them clearly** — e.g., `hg38_ensembl_115` or `grch37_refseq_110` — not `genome1`).

SSH into the server:

```bash
ssh root@YOUR_SERVER_IP
su - fastvep
mkdir -p /opt/fastvep/data
cd /opt/fastvep/data
```

Also upload the `fastvep` CLI binary (needed to build SA databases):

```bash
# From dev machine
scp target/x86_64-unknown-linux-gnu/release/fastvep root@YOUR_SERVER_IP:/opt/fastvep/bin/
```

Install `samtools` for FASTA indexing:

```bash
apt install -y samtools
```

---

### Human — GRCh38 / hg38 (Ensembl 115, primary assembly)

**This is the recommended default.** Most modern VCFs are aligned to GRCh38.

```bash
mkdir -p hg38_ensembl_115/sa && cd hg38_ensembl_115

# Gene annotation (~1.1GB compressed → ~3GB uncompressed)
wget https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.gff3.gz
gunzip Homo_sapiens.GRCh38.115.gff3.gz

# Reference FASTA (~850MB compressed → ~2.9GB uncompressed)
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa

cd /opt/fastvep/data
```

### Human — GRCh37 / hg19 (Ensembl 113, last GRCh37 release)

Still widely used for clinical and legacy data. Ensembl froze GRCh37 support at release 113.

```bash
mkdir -p hg19_ensembl_113/sa && cd hg19_ensembl_113

wget https://ftp.ensembl.org/pub/grch37/release-113/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz
gunzip Homo_sapiens.GRCh37.87.gff3.gz

wget https://ftp.ensembl.org/pub/grch37/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa

cd /opt/fastvep/data
```

### Human — CHM13 / T2T (GCA_009914755.4)

Telomere-to-telomere complete reference. Use if your VCFs are aligned to T2T-CHM13v2.0.

```bash
mkdir -p t2t_chm13v2/sa && cd t2t_chm13v2

# Gene annotation (CAT+Liftoff GFF3)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/vertebrate_mammalian/Homo_sapiens/T2T-CHM13v2.0/GCF_009914755.1-RS_2023_10/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz
gunzip GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz
mv GCF_009914755.1_T2T-CHM13v2.0_genomic.gff t2t_chm13v2.gff3

# Reference FASTA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz
gunzip GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz
mv GCA_009914755.4_T2T-CHM13v2.0_genomic.fna t2t_chm13v2.fa
samtools faidx t2t_chm13v2.fa

cd /opt/fastvep/data
```

### Human — GRCh38 (GENCODE 47 alternative to Ensembl)

If you prefer GENCODE's annotation over Ensembl's (they're nearly identical but numbered differently):

```bash
mkdir -p hg38_gencode_47/sa && cd hg38_gencode_47

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gff3.gz
gunzip gencode.v47.annotation.gff3.gz

# FASTA: reuse the Ensembl primary_assembly FASTA (same sequence), or download GENCODE's:
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
samtools faidx GRCh38.primary_assembly.genome.fa

cd /opt/fastvep/data
```

---

### Model organisms

All of these are available from Ensembl release 115. Same pattern: `gff3/<species>/` and `fasta/<species>/dna/`.

```bash
cd /opt/fastvep/data

# Mouse (Mus musculus) — GRCm39
mkdir -p mouse_grcm39_ensembl_115 && cd mouse_grcm39_ensembl_115
wget https://ftp.ensembl.org/pub/release-115/gff3/mus_musculus/Mus_musculus.GRCm39.115.gff3.gz
gunzip Mus_musculus.GRCm39.115.gff3.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
samtools faidx Mus_musculus.GRCm39.dna.primary_assembly.fa
cd /opt/fastvep/data

# Rat (Rattus norvegicus) — mRatBN7.2
mkdir -p rat_mratbn72_ensembl_115 && cd rat_mratbn72_ensembl_115
wget https://ftp.ensembl.org/pub/release-115/gff3/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.115.gff3.gz
gunzip Rattus_norvegicus.mRatBN7.2.115.gff3.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz
gunzip Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz
samtools faidx Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa
cd /opt/fastvep/data

# Zebrafish (Danio rerio) — GRCz11
mkdir -p zebrafish_grcz11_ensembl_115 && cd zebrafish_grcz11_ensembl_115
wget https://ftp.ensembl.org/pub/release-115/gff3/danio_rerio/Danio_rerio.GRCz11.115.gff3.gz
gunzip Danio_rerio.GRCz11.115.gff3.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
gunzip Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
samtools faidx Danio_rerio.GRCz11.dna.primary_assembly.fa
cd /opt/fastvep/data

# Fruit fly (Drosophila melanogaster) — BDGP6.46
mkdir -p drosophila_bdgp6_ensembl_115 && cd drosophila_bdgp6_ensembl_115
wget https://ftp.ensembl.org/pub/release-115/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.115.gff3.gz
gunzip Drosophila_melanogaster.BDGP6.46.115.gff3.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz
gunzip Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz
samtools faidx Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa
cd /opt/fastvep/data

# C. elegans — WBcel235
mkdir -p celegans_wbcel235_ensembl_115 && cd celegans_wbcel235_ensembl_115
wget https://ftp.ensembl.org/pub/release-115/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.115.gff3.gz
gunzip Caenorhabditis_elegans.WBcel235.115.gff3.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz
gunzip Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz
samtools faidx Caenorhabditis_elegans.WBcel235.dna.toplevel.fa
cd /opt/fastvep/data

# Yeast (S. cerevisiae) — R64-1-1
mkdir -p yeast_r64_ensembl_115 && cd yeast_r64_ensembl_115
wget https://ftp.ensembl.org/pub/release-115/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.115.gff3.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.115.gff3.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
samtools faidx Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
cd /opt/fastvep/data

# Arabidopsis — TAIR10
mkdir -p arabidopsis_tair10_ensembl_115 && cd arabidopsis_tair10_ensembl_115
wget https://ftp.ensemblgenomes.org/pub/plants/release-58/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.58.gff3.gz
gunzip Arabidopsis_thaliana.TAIR10.58.gff3.gz
wget https://ftp.ensemblgenomes.org/pub/plants/release-58/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
samtools faidx Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
cd /opt/fastvep/data
```

---

### Build supplementary annotation databases (human GRCh38 only)

SA databases are human-specific and build-specific. ClinVar/gnomAD/dbSNP have separate GRCh37 and GRCh38 VCFs — make sure you grab the right one for each genome subdirectory.

For **hg38_ensembl_115** (GRCh38):

```bash
cd /opt/fastvep/data/hg38_ensembl_115/sa

# ClinVar (~30MB download, ~40MB .osa) — updates monthly
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
/opt/fastvep/bin/fastvep sa-build --source clinvar \
  -i clinvar.vcf.gz -o clinvar --assembly GRCh38

# dbSNP (~20GB download, ~5GB .osa)
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
/opt/fastvep/bin/fastvep sa-build --source dbsnp \
  -i GCF_000001405.40.gz -o dbsnp --assembly GRCh38

# gnomAD v4 genomes (~25GB per chromosome, ~30GB .osa total)
for chr in {1..22} X Y; do
  wget "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz"
done
/opt/fastvep/bin/fastvep sa-build --source gnomad \
  -i gnomad.genomes.v4.1.sites.chr1.vcf.bgz -o gnomad --assembly GRCh38

# 1000 Genomes (~10GB download)
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
/opt/fastvep/bin/fastvep sa-build --source onekg \
  -i 1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
  -o onekg --assembly GRCh38

# REVEL (~2GB download)
# Manual download from https://sites.google.com/site/revelgenomics/downloads
/opt/fastvep/bin/fastvep sa-build --source revel \
  -i revel-v1.3_all_chromosomes.csv.zip -o revel --assembly GRCh38

# SpliceAI (~30GB download)
# Manual download from https://basespace.illumina.com/s/otSPW8hnhaZR
/opt/fastvep/bin/fastvep sa-build --source spliceai \
  -i spliceai_scores.raw.snv.hg38.vcf.gz -o spliceai --assembly GRCh38

# COSMIC (requires academic/commercial license)
# Register at https://cancer.sanger.ac.uk/cosmic/download

# Clean up downloaded source files
rm -f *.vcf.gz *.vcf.bgz *.csv.zip
```

For **hg19_ensembl_113** (GRCh37):

```bash
cd /opt/fastvep/data/hg19_ensembl_113/sa

# ClinVar GRCh37
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
/opt/fastvep/bin/fastvep sa-build --source clinvar \
  -i clinvar.vcf.gz -o clinvar --assembly GRCh37

# dbSNP GRCh37 (legacy — b151 is the last GRCh37 release)
wget https://ftp.ncbi.nih.gov/snp/pre_build152/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
/opt/fastvep/bin/fastvep sa-build --source dbsnp \
  -i 00-All.vcf.gz -o dbsnp --assembly GRCh37

# gnomAD v2.1.1 for GRCh37 (v3+ is GRCh38 only)
for chr in {1..22} X Y; do
  wget "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.${chr}.vcf.bgz"
done
/opt/fastvep/bin/fastvep sa-build --source gnomad \
  -i gnomad.genomes.r2.1.1.sites.1.vcf.bgz -o gnomad --assembly GRCh37

# (SpliceAI has separate hg19 and hg38 distributions)
# Manual download: spliceai_scores.raw.snv.hg19.vcf.gz
/opt/fastvep/bin/fastvep sa-build --source spliceai \
  -i spliceai_scores.raw.snv.hg19.vcf.gz -o spliceai --assembly GRCh37

rm -f *.vcf.gz *.vcf.bgz
```

---

### Disk space estimate per genome

Human GRCh38 with all SA databases:

| File | Size |
|------|------|
| GFF3 | 1.1 GB |
| GFF3 cache (auto-built on first run) | ~370 MB |
| FASTA + .fai | 2.9 GB |
| ClinVar .osa | ~40 MB |
| dbSNP .osa | ~5 GB |
| gnomAD .osa | ~30 GB |
| 1000 Genomes .osa | ~3 GB |
| REVEL .osa | ~1 GB |
| SpliceAI .osa | ~5 GB |
| **Total per human genome** | **~50 GB** |

Model organisms (without SA): typically 2-4 GB each. The smaller genomes (yeast, arabidopsis) are under 500 MB.

> **Disk planning**: For 2 human builds (hg38 + hg19) with full SA + a handful of model organisms, plan for ~110-130 GB. Use Hetzner CCX33 (160 GB SSD, ~$55/month) or attach a Hetzner Volume ($0.05/GB/month, e.g., 200GB = $10/month).

## 4. Configure the Service

### systemd unit

Two configurations depending on your setup:

**Option A: Single genome (simplest)**

Good if you only host one genome (e.g., just human GRCh38).

```bash
cat > /etc/systemd/system/fastvep.service << 'EOF'
[Unit]
Description=fastVEP Variant Annotation Server
After=network.target

[Service]
Type=simple
User=fastvep
Group=fastvep
WorkingDirectory=/opt/fastvep

Environment=FASTVEP_PORT=3000
Environment=FASTVEP_BIND=127.0.0.1
Environment=FASTVEP_GFF3=/opt/fastvep/data/hg38_ensembl_115/Homo_sapiens.GRCh38.115.gff3
Environment=FASTVEP_FASTA=/opt/fastvep/data/hg38_ensembl_115/Homo_sapiens.GRCh38.dna.primary_assembly.fa
Environment=FASTVEP_SA_DIR=/opt/fastvep/data/hg38_ensembl_115/sa

ExecStart=/opt/fastvep/bin/fastvep-web
Restart=always
RestartSec=5

LimitNOFILE=65535
MemoryMax=12G

[Install]
WantedBy=multi-user.target
EOF
```

**Option B: Multi-genome (recommended)**

Uses `--data-dir` so users can switch genomes from the web UI dropdown. Each subdirectory of `--data-dir` appears as a selectable genome.

```bash
cat > /etc/systemd/system/fastvep.service << 'EOF'
[Unit]
Description=fastVEP Variant Annotation Server
After=network.target

[Service]
Type=simple
User=fastvep
Group=fastvep
WorkingDirectory=/opt/fastvep

Environment=FASTVEP_PORT=3000
Environment=FASTVEP_BIND=127.0.0.1
# Data directory containing one subdirectory per genome (hg38_ensembl_115, hg19_ensembl_113, ...)
Environment=FASTVEP_DATA_DIR=/opt/fastvep/data
# The genome loaded at startup (users can switch via the web UI)
Environment=FASTVEP_GFF3=/opt/fastvep/data/hg38_ensembl_115/Homo_sapiens.GRCh38.115.gff3
Environment=FASTVEP_FASTA=/opt/fastvep/data/hg38_ensembl_115/Homo_sapiens.GRCh38.dna.primary_assembly.fa
Environment=FASTVEP_SA_DIR=/opt/fastvep/data/hg38_ensembl_115/sa

ExecStart=/opt/fastvep/bin/fastvep-web
Restart=always
RestartSec=5

LimitNOFILE=65535
# Allow more memory for multi-genome switching (each loaded genome uses ~4-6GB)
MemoryMax=24G

[Install]
WantedBy=multi-user.target
EOF
```

Start the service:

```bash
systemctl daemon-reload
systemctl enable fastvep
systemctl start fastvep

# Check it started
systemctl status fastvep
journalctl -u fastvep -f   # watch logs
```

First start takes ~10-15s (initial GFF3 parse + cache build). Subsequent starts are ~3-4s (from cache). Switching genomes in the UI takes ~3-4s.

### nginx reverse proxy

```bash
cat > /etc/nginx/sites-available/fastvep << 'EOF'
server {
    listen 80;
    server_name your-domain.com;

    # Rate limiting: 10 requests/sec per IP, burst 20
    limit_req_zone $binary_remote_addr zone=fastvep:10m rate=10r/s;

    location / {
        limit_req zone=fastvep burst=20 nodelay;

        proxy_pass http://127.0.0.1:3000;
        proxy_http_version 1.1;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;

        # Increase body size for large VCF uploads
        client_max_body_size 10M;

        # Timeouts for slow annotation requests
        proxy_read_timeout 120s;
        proxy_send_timeout 120s;
    }
}
EOF

ln -sf /etc/nginx/sites-available/fastvep /etc/nginx/sites-enabled/
rm -f /etc/nginx/sites-enabled/default
nginx -t && systemctl reload nginx
```

### TLS with Let's Encrypt (or Cloudflare)

**Option A: Let's Encrypt (free, auto-renewing)**

```bash
certbot --nginx -d your-domain.com
# Follow prompts; auto-renews via systemd timer
```

**Option B: Cloudflare (recommended — free DDoS protection + CDN)**

1. **Sign up** at [cloudflare.com](https://www.cloudflare.com/) (free plan)
2. **Add your domain**: Enter your domain, Cloudflare scans existing DNS records
3. **Update nameservers**: Point your domain registrar's nameservers to Cloudflare's (shown during setup)
4. **DNS records**: Add an A record:
   - Type: `A`
   - Name: `@` (or subdomain like `vep`)
   - IPv4: your Hetzner server IP
   - Proxy: **Proxied** (orange cloud icon ON) — this enables CDN + DDoS protection
5. **SSL/TLS settings** (Cloudflare dashboard → SSL/TLS):
   - Mode: **Full (strict)**
   - Edge Certificates: enabled (automatic)
6. **Origin certificate** (SSL/TLS → Origin Server):
   - Click "Create Certificate"
   - Let Cloudflare generate key + cert (RSA 2048, 15 years)
   - Save `origin.pem` and `origin-key.pem` to server:
   ```bash
   mkdir -p /etc/ssl/cloudflare
   # paste cert into /etc/ssl/cloudflare/origin.pem
   # paste key into /etc/ssl/cloudflare/origin-key.pem
   chmod 600 /etc/ssl/cloudflare/origin-key.pem
   ```
7. **Update nginx** to use the origin certificate:
   ```nginx
   server {
       listen 443 ssl;
       server_name your-domain.com;
       ssl_certificate /etc/ssl/cloudflare/origin.pem;
       ssl_certificate_key /etc/ssl/cloudflare/origin-key.pem;

       # ... same proxy_pass config as above ...
   }
   server {
       listen 80;
       server_name your-domain.com;
       return 301 https://$host$request_uri;
   }
   ```
8. **Recommended Cloudflare settings** (free tier):
   - Caching → Cache Rules: Cache the HTML page (`/` and `/index.html`) with short TTL (1 hour)
   - Security → WAF: Managed rules enabled (automatic)
   - Speed → Auto Minify: Enable HTML/CSS/JS minification
   - Analytics: Monitor traffic and threats for free

## 5. Verify

```bash
# From your machine
curl -s https://your-domain.com/api/status | python3 -m json.tool

# Should show:
# {
#     "status": "ok",
#     "transcripts": 508530,
#     "has_fasta": true,
#     "sa_sources": ["ClinVar", "dbSNP", "gnomAD", ...]
# }

# Test annotation
curl -s -X POST https://your-domain.com/api/annotate \
  -H 'Content-Type: application/json' \
  -d '{"vcf": "1\t65565\t.\tA\tG\t50\tPASS\t.", "pick": false}' | python3 -m json.tool
```

## 6. Maintenance

### Updates

```bash
# Build new release on dev machine
cross build --release --target x86_64-unknown-linux-gnu -p fastvep-web

# Deploy (< 1 minute downtime)
scp target/x86_64-unknown-linux-gnu/release/fastvep-web root@YOUR_SERVER_IP:/opt/fastvep/bin/fastvep-web.new
ssh root@YOUR_SERVER_IP '
  mv /opt/fastvep/bin/fastvep-web.new /opt/fastvep/bin/fastvep-web
  chmod +x /opt/fastvep/bin/fastvep-web
  systemctl restart fastvep
'
```

### Update SA databases

```bash
# ClinVar updates monthly
ssh root@YOUR_SERVER_IP
cd /tmp
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
/opt/fastvep/bin/fastvep sa-build --source clinvar \
  -i clinvar.vcf.gz -o /opt/fastvep/data/sa/clinvar
systemctl restart fastvep
```

### Monitoring

```bash
# Check service health
systemctl status fastvep

# Watch logs
journalctl -u fastvep -f

# Memory usage
ps aux | grep fastvep-web

# Disk usage
du -sh /opt/fastvep/data/*
```

### Backups

The only state that matters is the SA databases (everything else is re-downloadable):

```bash
# Weekly backup of SA databases (add to crontab)
tar czf /opt/fastvep/backup/sa_$(date +%Y%m%d).tar.gz /opt/fastvep/data/sa/
```

## 7. Scaling (if needed later)

| Situation | Solution | Cost |
|-----------|----------|------|
| More traffic | Upgrade to CCX33 (8 vCPU, 32GB) | ~$55/month |
| High availability | Second server + Cloudflare load balancing | ~$60/month |
| Multiple genomes | Use `--data-dir` with subdirectories per organism | Same server |
| Batch WGS jobs | Run `fastvep` CLI in background, rate-limited | Same server |
| Global users | Add a second server in another region | ~$60/month |

## Architecture Diagram

```
Users (browser)
      │
      ▼
  Cloudflare (free)
  ├─ DDoS protection
  ├─ TLS termination
  └─ Caching
      │
      ▼
  Hetzner CCX23 ($30/mo)
  ├─ nginx (reverse proxy, rate limiting)
  └─ fastvep-web (port 3000)
     ├─ 508k transcripts (GRCh38, Ensembl 115)
     ├─ FASTA reference (memory-mapped, 2.9GB)
     └─ SA databases
        ├─ ClinVar
        ├─ gnomAD v4
        ├─ dbSNP
        ├─ 1000 Genomes
        ├─ REVEL
        └─ SpliceAI
```
