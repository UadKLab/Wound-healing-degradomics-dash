# Wound Healing Degradomics

Dash Python dashboard.

Look at `.env-example` to set up environment variables in an local `.env` file. In order to run the production code you will need a S3 bucket to store the same data as in `data` folder, along with the necessary authentication to access the S3 bucket.

Run instructions for development:   
// for creating
`docker compose -f compose.dev.yaml up`
// for stopping
`docker compose -f compose.dev.yaml down`
