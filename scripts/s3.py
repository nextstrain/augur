"""
Script to sync files between local disk and S3 buckets.

# Go to flu build.
cd augur/builds/flu

# Download flu H3N2 data into auspice directory.
python ../../scripts/s3.py pull dev-data \
    --prefixes flu_h3n2 --to auspice

# Upload flu H3N2 data to S3 dev bucket.
python ../../scripts/s3.py push dev-data \
    auspice/flu_h3n2_*

# Sync H3N2 data from one bucket to another and create CloudFront invalidation.
python ../../scripts/s3.py sync dev-data production-data \
    --prefixes flu_h3n2
"""
import argparse
import boto3
import botocore
import gzip
import io
import logging
import os
import shutil
import time

# Map S3 buckets to their corresponding CloudFront ids.
CLOUDFRONT_ID_BY_BUCKET = {
    "nextstrain-staging": "E3L83FTHWUN0BV",
    "nextstrain-data": "E3LB0EWZKCCV"
}


def get_bucket_keys_by_prefixes(bucket, prefixes):
    """Get objects from the given bucket instance that match the given list of
    prefixes. If no prefixes are given, get all objects.

    Args:
        bucket: S3 Bucket instance
        prefixes: NoneType or list of key prefixes to filter objects in the bucket by

    Returns:
        list: sorted list of object keys matching the given prefixes
    """
    if prefixes is not None:
        object_keys = []
        for prefix in prefixes:
            keys = [obj.key for obj in bucket.objects.filter(Prefix=prefix)]
            object_keys.extend(keys)

        object_keys = sorted(set(object_keys))
    else:
        object_keys = sorted([obj.key for obj in bucket.objects.all()])

    return object_keys


def create_cloudfront_invalidation(bucket_name, keys):
    """Create a cache invalidation for the given files if a CloudFront id is given.

    Args:
        bucket_name: an S3 bucket name that may or may not have a CloudFront id
        keys: a list of keys to invalidate cache for

    Returns:
        dict or NoneType: CloudFront API response or None if no CloudFront is defined for the given bucket
    """
    # Setup logging.
    logger = logging.getLogger(__name__)

    # Find CloudFront id for the given bucket name.
    cloudfront_id = CLOUDFRONT_ID_BY_BUCKET.get(bucket_name)
    if cloudfront_id is None:
        logger.warning("Could not find a CloudFront id for the S3 bucket '%s'" % bucket_name)
        return

    logger.info("Creating invalidation for %i keys in the CloudFront distribution '%s'" % (len(keys), cloudfront_id))

    # Connect to CloudFront.
    cloudfront = boto3.client("cloudfront")

    # Create the invalidation batch. Top-level keys require a "/" prefix for
    # proper invalidation.
    invalidation_batch = {
        "Paths": {
            "Quantity": len(keys),
            "Items": ["/%s" % key for key in keys]
        },
        "CallerReference": str(time.time())
    }
    logger.debug("Invalidation batch: %s" % str(invalidation_batch))

    # Create the invalidation.
    response = cloudfront.create_invalidation(
        DistributionId=cloudfront_id,
        InvalidationBatch=invalidation_batch
    )
    logger.debug("CloudFront response: %s" % str(response))

    return response


def push(bucket_name, files, dryrun=False):
    """Push the given files to the given S3 bucket and optionally invalidate the
    cache for a given CloudFront id.

    Args:
        bucket_name: S3 bucket to pull from
        files: a list of local files to upload to the given bucket
        dryrun: boolean indicating whether files should be downloaded or not
    """
    # Setup logging.
    logger = logging.getLogger(__name__)

    # Create a distinct list of files to push.
    files = list(set(files))

    # Confirm that all given file paths are proper files.
    for file_name in files:
        assert os.path.isfile(file_name), "The requested input '%s' is not a proper file" % file_name

    # Connect to S3.
    s3 = boto3.resource("s3")
    bucket = s3.Bucket(bucket_name)

    # Upload local files, stripping directory names from the given file paths
    # for the S3 keys.
    logger.info("Uploading %i files to bucket '%s'" % (len(files), bucket_name))
    s3_keys = []
    for file_name in files:
        s3_key = os.path.split(file_name)[-1]
        s3_keys.append(s3_key)
        logger.debug("Uploading '%s' as '%s'" % (file_name, s3_key))

        if not dryrun:
            # Open uncompressed file to be uploaded.
            with open(file_name, "rb") as fh:
                # Compress the file in memory.
                compressed_fh = io.BytesIO()
                with gzip.GzipFile(fileobj=compressed_fh, mode="wb") as gz:
                    shutil.copyfileobj(fh, gz)

                compressed_fh.seek(0)

                # Upload the compressed file with associated metadata.
                bucket.upload_fileobj(
                    compressed_fh,
                    s3_key,
                    {"ContentType": "gzip", "ContentEncoding": "application/json"}
                )

    # Create a CloudFront invalidation for the destination bucket.
    if not dryrun:
        response = create_cloudfront_invalidation(bucket_name, s3_keys)


def pull(bucket_name, prefixes=None, local_dir=None, dryrun=False):
    """Pull files from the given S3 bucket. Optionally, only pull files that match
    the given list of filename prefixes.

    Args:
        bucket_name: S3 bucket to pull from
        prefixes: a list of key prefixes to filter objects in the bucket by
        local_dir: a local directory to download files into
        dryrun: boolean indicating whether files should be downloaded or not
    """
    # Setup logging.
    logger = logging.getLogger(__name__)

    # Confirm that the given local directory is a real directory.
    if local_dir is not None:
        assert os.path.isdir(local_dir), "The requested output directory '%s' is not a proper directory." % local_dir

    # Connect to S3.
    s3 = boto3.resource("s3")

    # Get a list of all objects in the requested bucket.
    bucket = s3.Bucket(bucket_name)

    # Get keys by prefixes.
    object_keys = get_bucket_keys_by_prefixes(bucket, prefixes)

    # Download objects.
    logger.info("Downloading %i files from bucket '%s'" % (len(object_keys), bucket_name))
    for key in object_keys:
        # Download into a local directory if requested.
        if local_dir is not None:
            local_key = os.path.join(local_dir, key)
        else:
            local_key = key

        logger.debug("Downloading '%s' as '%s'" % (key, local_key))
        if not dryrun:
            with open(local_key, "wb") as fh:
                # Download the compressed file into memory.
                compressed_fh = io.BytesIO()
                bucket.download_fileobj(key, compressed_fh)
                compressed_fh.seek(0)

                # Write the uncompressed data from the file to disk.
                with gzip.GzipFile(fileobj=compressed_fh, mode="rb") as gz:
                    shutil.copyfileobj(gz, fh)


def sync(source_bucket_name, destination_bucket_name, prefixes=None, dryrun=False):
    """Sync files from a given source bucket to a given destination bucket. An
    optional list of prefixes will restrict the files being synced between
    buckets.

    Args:
        source_bucket_name: name of bucket to copy files from
        destination_bucket_name: name of bucket to copy files to
        prefixes: a list of key prefixes to filter objects in the bucket by
        dryrun: boolean indicating whether files should be downloaded or not
    """
    # Setup logging.
    logger = logging.getLogger(__name__)

    # Connect to S3.
    s3 = boto3.resource("s3")

    # Get source bucket.
    source_bucket = s3.Bucket(source_bucket_name)

    # Get source bucket keys by prefixes.
    object_keys = get_bucket_keys_by_prefixes(source_bucket, prefixes)
    logger.info("Syncing %i files from '%s' to '%s'" % (len(object_keys), source_bucket_name, destination_bucket_name))

    # Copy objects from source to destination by key.
    for key in object_keys:
        logger.debug("Copying '%s'" % key)

        if not dryrun:
            copy_source = {
                'Bucket': source_bucket_name,
                'Key': key
            }
            s3.meta.client.copy(copy_source, destination_bucket_name, key)

    # Create a CloudFront invalidation for the destination bucket.
    if not dryrun:
        response = create_cloudfront_invalidation(destination_bucket_name, object_keys)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", "-v", action="store_const", dest="loglevel", const=logging.INFO, help="Enable verbose logging")
    parser.add_argument("--debug", "-d", action="store_const", dest="loglevel", const=logging.DEBUG, help="Enable debugging logging")
    parser.add_argument("--dryrun", "-n", action="store_true", help="Perform a dryrun without uploading or downloading any files")

    subparsers = parser.add_subparsers(dest="command_name")

    parser_push = subparsers.add_parser("push")
    parser_push.add_argument("bucket", help="S3 bucket to push files to")
    parser_push.add_argument("files", nargs="+", help="One or more sets of files to push to the given bucket")
    parser_push.set_defaults(func=push)

    parser_pull = subparsers.add_parser("pull")
    parser_pull.add_argument("bucket", help="S3 bucket to pull files from")
    parser_pull.add_argument("--prefixes", "-p", nargs="+", help="One or more file prefixes to match in the given bucket")
    parser_pull.add_argument("--local_dir", "--to", "-t", help="Local directory to download files into")
    parser_pull.set_defaults(func=pull)

    parser_sync = subparsers.add_parser("sync")
    parser_sync.add_argument("source_bucket", help="Source S3 bucket")
    parser_sync.add_argument("destination_bucket", help="Destination S3 bucket")
    parser_sync.add_argument("--prefixes", "-p", nargs="+", help="One or more prefixes for files to sync between buckets")
    parser_sync.set_defaults(func=sync)

    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel)

    try:
        if args.command_name == "push":
            args.func(args.bucket, args.files, args.dryrun)
        elif args.command_name == "pull":
            args.func(args.bucket, args.prefixes, args.local_dir, args.dryrun)
        elif args.command_name == "sync":
            args.func(args.source_bucket, args.destination_bucket, args.prefixes, args.dryrun)
    except botocore.exceptions.NoCredentialsError, e:
        parser.error("Unable to locate AWS credentials. Set environment variables for AWS_SECRET_ACCESS_KEY and AWS_ACCESS_KEY_ID.")
    except Exception, e:
        parser.error(e.message)
