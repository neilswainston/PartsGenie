#!/usr/bin/bash

# Set constants:
PROJECT_ID=personal-parts-genie
SERVICE_ACCOUNT_NAME=${PROJECT_ID}-sa
LOCATION=europe-west1
ARTIFACT_REPO_NAME=${PROJECT_ID}-repo
ARTIFACT_CLEANUP_POLICY=cleanup-policy.json
TRIGGER_NAME=${PROJECT_ID}-build-trigger
GITHUB_REPO_NAME=PartsGenie
GITHUB_REPO_OWNER=neilswainston
BUILD_CONFIG=cloudbuild.yaml
BRANCH_PATTERN="^main$"

# Enable services:
gcloud services enable \
    --project ${PROJECT_ID} \
    artifactregistry.googleapis.com \
    cloudbuild.googleapis.com \
    run.googleapis.com

# Create service account:
gcloud iam service-accounts create ${SERVICE_ACCOUNT_NAME} \
    --project ${PROJECT_ID}

# Add roles to service account:
gcloud projects add-iam-policy-binding ${PROJECT_ID} \
    --member="serviceAccount:${SERVICE_ACCOUNT_NAME}@${PROJECT_ID}.iam.gserviceaccount.com" \
    --role="roles/iam.serviceAccountUser"

gcloud projects add-iam-policy-binding ${PROJECT_ID} \
    --member="serviceAccount:${SERVICE_ACCOUNT_NAME}@${PROJECT_ID}.iam.gserviceaccount.com" \
    --role="roles/artifactregistry.writer"

gcloud projects add-iam-policy-binding ${PROJECT_ID} \
    --member="serviceAccount:${SERVICE_ACCOUNT_NAME}@${PROJECT_ID}.iam.gserviceaccount.com" \
    --role="roles/logging.logWriter"

gcloud projects add-iam-policy-binding ${PROJECT_ID} \
    --member="serviceAccount:${SERVICE_ACCOUNT_NAME}@${PROJECT_ID}.iam.gserviceaccount.com" \
    --role="roles/run.developer"

gcloud projects add-iam-policy-binding ${PROJECT_ID} \
    --member="serviceAccount:${SERVICE_ACCOUNT_NAME}@${PROJECT_ID}.iam.gserviceaccount.com" \
    --role="roles/storage.admin"

# Create Artefact Registry for Docker images:
gcloud artifacts repositories create ${ARTIFACT_REPO_NAME} \
    --project ${PROJECT_ID} \
    --location ${LOCATION} \
    --repository-format docker

gcloud artifacts repositories set-cleanup-policies ${ARTIFACT_REPO_NAME} \
    --project ${PROJECT_ID} \
    --location ${LOCATION} \
    --policy ${ARTIFACT_CLEANUP_POLICY} \
    --no-dry-run

# Set up Cloud Build:
gcloud builds triggers create github \
    --name ${TRIGGER_NAME} \
    --project ${PROJECT_ID} \
    --region ${LOCATION} \
    --repo-name ${GITHUB_REPO_NAME} \
    --repo-owner ${GITHUB_REPO_OWNER} \
    --build-config ${BUILD_CONFIG} \
    --branch-pattern ${BRANCH_PATTERN} \
    --service-account projects/${PROJECT_ID}/serviceAccounts/${SERVICE_ACCOUNT_NAME}@${PROJECT_ID}.iam.gserviceaccount.com