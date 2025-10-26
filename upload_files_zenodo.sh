#!/bin/bash

ZENODO_TOKEN=SVnlsedTXNBbaGQKgWJmXmI6vF0ONzaj295F7G7aZyGQyFaqW1hI8EzrKLL9


API_BASE="https://zenodo.org/api"

# 1) Create a new draft deposition
DEPOSITION_JSON=$(curl -fsS \
  -H "Authorization: Bearer $ZENODO_TOKEN" \
  -H "Content-Type: application/json" \
  -X POST "$API_BASE/deposit/depositions" \
  -d '{}')

DEPOSITION_ID=$(echo "$DEPOSITION_JSON" | jq -r '.id')
BUCKET_URL=$(echo "$DEPOSITION_JSON" | jq -r '.links.bucket')

echo "Created deposition: ID=$DEPOSITION_ID"
echo "Bucket: $BUCKET_URL"

# 2) Upload all .sif files from envs/ (send only the filename to the bucket URL)
shopt -s nullglob
for FILE in envs/*.sif; do
  NAME=$(basename "$FILE")
  echo "Uploading $FILE as $NAME ..."
  curl -fS --progress-bar \
    -H "Authorization: Bearer $ZENODO_TOKEN" \
    --upload-file "$FILE" \
    "$BUCKET_URL/$NAME"
  echo
done

echo "All uploads done."

# 3) Upload the pbmc10x.h5mu dataset file
DATASET_FILE="data/datasets/pbmc10x/pbmc10x.h5mu"
DATASET_NAME=$(basename "$DATASET_FILE")
echo "Uploading dataset file $DATASET_FILE as $DATASET_NAME ..."
curl -fS --progress-bar \
  -H "Authorization:  Bearer $ZENODO_TOKEN" \
  --upload-file "$DATASET_FILE" \
  "$BUCKET_URL/$DATASET_NAME"
echo
echo "Dataset file upload done."


# 4) Update the deposition metadata (copy from an existing deposition)
curl -s \
  -H "Authorization: Bearer $ZENODO_TOKEN" \
  https://zenodo.org/api/deposit/depositions/17444910\
  | jq '.metadata' > metadata.json

curl -s \
  -H "Authorization: Bearer $ZENODO_TOKEN" \
  -H "Content-Type: application/json" \
  -X PUT https://zenodo.org/api/deposit/depositions/$DEPOSITION_ID \
  -d "{\"metadata\": $(cat metadata.json)}"