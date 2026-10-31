#!/bin/bash
# Copy a built site to Stanford personal web space (AFS ~/WWW, served by web.stanford.edu),
# readable only after Stanford's web login: an .htaccess file limits it to people with a
# Stanford affiliation (UIT's "Stanford only" WebAuth block, https://uit.stanford.edu/node/2272).
# The pages hold the unpublished proteome's names and hits and the 51 drawn proteins'
# sequences, so the site must not be public.
#
#   protein-pages/deploy_stanford.sh <site folder> [SUNet ID] [folder name]
#
# Uploads through rice.stanford.edu (FarmShare), the machine that mounts AFS; it asks for the
# SUNet password and Duo. The site then opens at https://web.stanford.edu/~<SUNet ID>/<folder>/.
set -euo pipefail
SITE=${1%/}
ID=${2:-olgabot}
DIR=${3:-botryllus-proteins}
[ -f "$SITE/index.html" ] || { echo "no index.html in $SITE" >&2; exit 1; }
cat > "$SITE/.htaccess" <<'HT'
AuthType WebAuth
require privgroup stanford:stanford
require privgroup stanford:administrative
require privgroup stanford:academic
HT
# aklog turns the login's Kerberos ticket into the AFS token rsync needs to write to ~/WWW
# macOS's own rsync (openrsync) fails against Linux rsync servers, so this uses pixi's rsync 3
pixi exec rsync -a --delete --exclude .DS_Store --stats --rsync-path="aklog && mkdir -p afs-home/WWW/$DIR && rsync" \
  "$SITE/" "$ID@rice.stanford.edu:afs-home/WWW/$DIR/"
echo "https://web.stanford.edu/~$ID/$DIR/"
