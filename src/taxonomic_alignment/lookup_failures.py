"""Tell a failed WoRMS or GBIF request apart from a name the backbone does not contain."""

import logging
import random
import time

# One try plus two retries. Enough for a short outage or a rate limit, without
# turning a dead connection into a very long run.
LOOKUP_MAX_ATTEMPTS = 3


class TaxonomicLookupFailed(RuntimeError):
    """
    A taxonomic API request still failed after retries.
    Callers must not record these names as incertae sedis.
    """


def call_with_retry(func, description):
    """
    Run an API call and retry it when it raises.

    A timeout, a rate limit, or a dropped connection raises. A name the
    backbone does not contain comes back as an empty result and does not raise,
    so it is not retried here.
    """
    last_error = None
    for attempt in range(1, LOOKUP_MAX_ATTEMPTS + 1):
        try:
            return func()
        except Exception as exc:
            last_error = exc
            logging.warning(
                "%s failed (attempt %s/%s): %s",
                description,
                attempt,
                LOOKUP_MAX_ATTEMPTS,
                exc,
            )
            if attempt < LOOKUP_MAX_ATTEMPTS:
                # Back off, plus a little jitter so parallel workers do not
                # all hit the API again at the same moment.
                time.sleep((2 ** attempt) + random.uniform(0, 0.5))
    raise last_error


def format_lookup_failure(api_name, failed_names):
    """One message for the log, the HTML report, and the raised error."""
    unique = []
    seen = set()
    for name in failed_names:
        text = "" if name is None else str(name).strip()
        if not text or text in seen:
            continue
        seen.add(text)
        unique.append(text)

    examples = ", ".join(unique[:5])
    if len(unique) > 5:
        examples = f"{examples} (+{len(unique) - 5} more)"

    return (
        f"{api_name} request failed for {len(unique)} taxon name(s) after {LOOKUP_MAX_ATTEMPTS} attempts. "
        "These names were not assigned to incertae sedis, because a failed request is not the same as a name the backbone does not contain. "
        "Failed requests were not cached, so the next run will try them again. "
        f"Examples: {examples}"
    )


def higher_classification_warning(api_name, count):
    """Warning text when the name matched but the lineage request did not."""
    count = int(count or 0)
    if count <= 0:
        return None
    return (
        f"{api_name} higherClassification lookup failed for {count} taxon(s) after {LOOKUP_MAX_ATTEMPTS} attempts. "
        "The scientific name was still assigned. The missing lineage was not cached, so the next run will try again."
    )


def empty_gbif_cache_warning(count):
    """Warning text when old empty GBIF cache entries are looked up again."""
    count = int(count or 0)
    if count <= 0:
        return None
    return (
        f"Ignored {count} empty entries in the previous GBIF match cache. "
        "Failed requests used to be stored as if the name did not match, so those names will be looked up again."
    )
