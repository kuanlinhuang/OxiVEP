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

SSH into the server:

```bash
ssh root@YOUR_SERVER_IP
su - fastvep
cd /opt/fastvep/data
```

### Reference genome and gene annotation (~4GB)

```bash
# GFF3 gene models (Ensembl 115)
wget https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.gff3.gz
gunzip Homo_sapiens.GRCh38.115.gff3.gz

# Reference FASTA
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Create FASTA index (required for memory-mapped access)
apt install -y samtools   # or download from htslib
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

### Build supplementary annotation databases

Upload the `fastvep` CLI binary too (or build on server):

```bash
# Upload CLI binary (on dev machine)
scp target/x86_64-unknown-linux-gnu/release/fastvep root@YOUR_SERVER_IP:/opt/fastvep/bin/
```

Back on the server:

```bash
mkdir -p /opt/fastvep/data/sa
cd /tmp

# --- ClinVar (~30MB download, ~40MB .osa) ---
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
/opt/fastvep/bin/fastvep sa-build --source clinvar \
  -i clinvar.vcf.gz -o /opt/fastvep/data/sa/clinvar

# --- dbSNP (~20GB download, ~5GB .osa) ---
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
/opt/fastvep/bin/fastvep sa-build --source dbsnp \
  -i GCF_000001405.40.gz -o /opt/fastvep/data/sa/dbsnp

# --- gnomAD v4 genomes (~25GB per chromosome, ~30GB .osa total) ---
# Download the joint sites VCF from https://gnomad.broadinstitute.org/downloads#v4
# Or process per-chromosome:
for chr in {1..22} X Y; do
  wget "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz"
done
# Concatenate and build (or build per-chromosome and merge):
/opt/fastvep/bin/fastvep sa-build --source gnomad \
  -i gnomad.genomes.v4.1.sites.chr1.vcf.bgz -o /opt/fastvep/data/sa/gnomad

# --- 1000 Genomes (~10GB download) ---
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
/opt/fastvep/bin/fastvep sa-build --source onekg \
  -i 1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
  -o /opt/fastvep/data/sa/onekg

# --- REVEL (~2GB download) ---
# Download from https://sites.google.com/site/revelgenomics/downloads
/opt/fastvep/bin/fastvep sa-build --source revel \
  -i revel-v1.3_all_chromosomes.csv.zip -o /opt/fastvep/data/sa/revel

# --- SpliceAI (~30GB download) ---
# Download from https://basespace.illumina.com/s/otSPW8hnhaZR
/opt/fastvep/bin/fastvep sa-build --source spliceai \
  -i spliceai_scores.raw.snv.hg38.vcf.gz -o /opt/fastvep/data/sa/spliceai

# --- COSMIC (requires license) ---
# Register at https://cancer.sanger.ac.uk/cosmic/download
# Download CosmicCodingMuts.vcf.gz
# /opt/fastvep/bin/fastvep sa-build --source cosmic \
#   -i CosmicCodingMuts.vcf.gz -o /opt/fastvep/data/sa/cosmic

# Clean up downloaded source files to reclaim disk
rm -f /tmp/*.vcf.gz /tmp/*.vcf.bgz
```

### Disk space estimate

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
| **Total on disk** | **~50 GB** |

> The 80GB SSD on CCX23 is sufficient. If you add all databases AND keep the source VCFs around, upgrade to CCX33 (160GB) or attach a volume.

## 4. Configure the Service

### systemd unit

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
Environment=FASTVEP_GFF3=/opt/fastvep/data/Homo_sapiens.GRCh38.115.gff3
Environment=FASTVEP_FASTA=/opt/fastvep/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa
Environment=FASTVEP_SA_DIR=/opt/fastvep/data/sa

ExecStart=/opt/fastvep/bin/fastvep-web
Restart=always
RestartSec=5

# Resource limits
LimitNOFILE=65535
MemoryMax=12G

[Install]
WantedBy=multi-user.target
EOF

systemctl daemon-reload
systemctl enable fastvep
systemctl start fastvep

# Check it started
systemctl status fastvep
journalctl -u fastvep -f   # watch logs
```

First start takes ~4 seconds (builds transcript cache). Subsequent starts load from cache.

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

**Option B: Cloudflare (free, DDoS protection)**

1. Add your domain to Cloudflare (free plan)
2. Point DNS A record to your server IP (proxy enabled, orange cloud)
3. SSL/TLS → Full (strict) in Cloudflare dashboard
4. Origin certificate: Create one in Cloudflare and install in nginx

Cloudflare is recommended — free DDoS protection, caching of static assets, and analytics.

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
