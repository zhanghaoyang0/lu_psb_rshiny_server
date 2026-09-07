# Builds zhanghaoyang0/rshiny-pon_del, a superset of zhanghaoyang0/rshiny
# with transvar installed directly (no sibling container / docker.sock needed).
FROM zhanghaoyang0/rshiny:latest

RUN pip3 install --no-cache-dir transvar==2.5.9.20190830 future

# transvar's own annotation database (hg38.*.transvardb etc.), baked in from
# the upstream transvar image so it doesn't need to be supplied externally.
COPY --from=zhouwanding/transvar:latest /anno /anno

ENV TRANSVAR_CFG=/anno/transvar.cfg
ENV TRANSVAR_DOWNLOAD_DIR=/anno/
