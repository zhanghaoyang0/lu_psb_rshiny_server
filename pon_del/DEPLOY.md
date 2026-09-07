# Deploying pon_del on a new machine

The `zhanghaoyang0/rshiny-pon_del` image is a superset of `zhanghaoyang0/rshiny`
with `transvar` installed directly (no sibling container, no
`/var/run/docker.sock` mount). It bakes in transvar's own annotation
database, but it does **not** bake in the app code, the protein DB, or the
reference genome — those are supplied via bind mounts at run time, same as
before.

## 1. Get the image

```bash
docker pull zhanghaoyang0/rshiny-pon_del:latest
```

Or build it locally from `pon_del/Dockerfile`:

```bash
docker build -t zhanghaoyang0/rshiny-pon_del:latest -f pon_del/Dockerfile pon_del/
```

## 2. Get the three things the image does NOT contain

| What | Size | Source | Mounted to |
|---|---|---|---|
| App code (this repo) | ~5.3G (all apps) | `git clone git@github.com:zhanghaoyang0/lu_psb_rshiny_server.git` | `/srv/shiny-server` |
| Protein DB | ~22G | copy from old host's `/home/ha0214zh/db/protein` (rsync/scp) | `/srv/shiny-server/db` |
| Reference genome (hg38.fa/hg19.fa + BLAST indices) | ~6.8G | copy from old host's `/home/ha0214zh/soft/vep/genome` | `/ref` |

## 3. Run

```bash
docker run -itd -p 8503:3838 \
  --user root \
  --name rshiny \
  -v <path-to-repo>:/srv/shiny-server \
  -v <path-to-protein-db>:/srv/shiny-server/db \
  -v <path-to-genome-dir>:/ref \
  zhanghaoyang0/rshiny-pon_del:latest
```

## 4. Verify

```bash
curl -s -o /dev/null -w '%{http_code}\n' http://localhost:8503/pon_del/   # expect 200
docker exec rshiny bash -c "transvar panno -i 'ENST00000242592.9:p.104_104del' --ensembl --reference /ref/hg38.fa"
```
