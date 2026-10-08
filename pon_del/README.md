# PON-Del has moved

The PON-Del web app is maintained in the private `zhanghaoyang0/pon_del` repository,
under `code/deploy` (deployment steps in its `DEPLOY.md`).

On the server this directory is only the mount point: the `rshiny` container mounts
`pon_del/code/deploy` over `/srv/shiny-server/pon_del`. Keep the directory, do not delete it.
