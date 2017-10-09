"""
Script to sync local files to S3 and between S3 buckets.
"""
import argparse
import boto3
import logging
import os

# Setup logger and handlers.
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger_handler = logging.StreamHandler()
logger_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger_handler.setFormatter(formatter)
logger.addHandler(logger_handler)


def push(bucket_name, files, cloudfront_id=None):
    """Push the given files to the given S3 bucket and optionally invalidate the
    cache for a given CloudFront id.
    """
    # Create a distinct list of files to push.
    files = list(set(files))

    # Connect to S3.
    s3 = boto3.resource("s3")
    bucket = s3.Bucket(bucket_name)

    # Upload local files, stripping directory names from the given file paths
    # for the S3 keys.
    logger.info("Uploading %i files to bucket '%s'" % (len(files), bucket_name))
    for file_name in files:
        s3_key = os.path.split(file_name)[-1]
        logger.debug("Uploading '%s' as '%s'" % (file_name, s3_key))
        bucket.upload_file(file_name, s3_key)

    if cloudfront_id is not None:
        logger.debug("Invalidating cache for CloudFront id '%s'" % cloudfront_id)


def pull(bucket_name, prefixes=None):
    """Pull files from the given S3 bucket. Optionally, only pull files that match
    the given list of filename prefixes.
    """
    # Connect to S3.
    s3 = boto3.resource("s3")

    # Get a list of all objects in the requested bucket.
    bucket = s3.Bucket(bucket_name)

    # Get objects that match the given list of prefixes. If no prefixes are
    # given, get all objects.
    if prefixes is not None:
        object_keys = []
        for prefix in prefixes:
            keys = [obj.key for obj in bucket.objects.filter(Prefix=prefix)]
            object_keys.extend(keys)

        object_keys = sorted(set(object_keys))
    else:
        object_keys = sorted([obj.key for obj in bucket.objects.all()])

    # Download objects.
    logger.info("Downloading %i objects" % (len(object_keys),))
    for key in object_keys:
        logger.debug("Downloading %s" % (key,))
        bucket.download_file(key, key)


def sync(source_bucket, destination_bucket, files, cloudfront_id=None):
    """Sync the given files from the given source S3 bucket to the given destination
    bucket and optionally invalidate the cache for a given CloudFront id.
    """
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")

    subparsers = parser.add_subparsers(dest="command_name")

    parser_push = subparsers.add_parser("push")
    parser_push.add_argument("--cloudfront_id", "-c", help="CloudFront id to use to create a cache invalidation")
    parser_push.add_argument("bucket", help="S3 bucket to push files to")
    parser_push.add_argument("files", nargs="+", help="One or more sets of files to push to the given bucket")
    parser_push.set_defaults(func=push)

    parser_pull = subparsers.add_parser("pull")
    parser_pull.add_argument("bucket", help="S3 bucket to pull files from")
    parser_pull.add_argument("--prefixes", nargs="+", help="One or more file prefixes to match in the given bucket")
    parser_pull.set_defaults(func=pull)

    parser_sync = subparsers.add_parser("sync")
    parser_sync.add_argument("--cloudfront_id", "-c", help="CloudFront id to use to create a cache invalidation")
    parser_sync.add_argument("source_bucket", help="S3 bucket to sync files from")
    parser_sync.add_argument("destination_bucket", help="S3 bucket to sync files to")
    parser_sync.add_argument("files", nargs="+", help="One or more sets of files to sync between the given buckets")
    parser_sync.set_defaults(func=sync)

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)
        logger_handler.setLevel(logging.DEBUG)

    if args.command_name == "push":
        args.func(args.bucket, args.files, args.cloudfront_id)
    elif args.command_name == "pull":
        args.func(args.bucket, args.prefixes)
    elif args.command_name == "sync":
        args.func(args.source_bucket, args.destination_bucket, args.files, args.cloudfront_id)
    else:
        args.func(args)
