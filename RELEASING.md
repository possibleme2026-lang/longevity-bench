# Releasing Longevity-Bench

> **Inherited document — upstream org details below are being adapted.**
> This release process comes from
> [Terminal-Bench-Science](https://github.com/harbor-framework/terminal-bench-science)
> along with the rest of the CI. The mechanics — semantic-versioned GitHub
> releases, `CITATION.cff` as the metadata source, the Zenodo/Software Heritage
> archival flow, the `tools/release-notes/generate.py` draft — are the same
> machinery this repository runs, so the steps and gates are accurate here.
> What is *not* accurate yet is anything naming an organisation or an external
> record: the `harbor-framework` org, the repo Zenodo is asked for access to, and
> the Zenodo concept record ID are upstream's. They are marked with
> `TODO(owner)` below. **This repository also has no `CITATION.cff` and no Zenodo
> deposit yet**, so steps 3–4 below and `.github/workflows/update-citation.yml`
> cannot run as written — create the citation file and your own Zenodo record
> first, or drop the Zenodo steps entirely.

Longevity-Bench releases are immutable benchmark snapshots. Each release
is published as a semantic-versioned GitHub release and, once the repository is
enabled in Zenodo, archived automatically with a DOI.

## Versioning

Release tags use `vMAJOR.MINOR.PATCH`:

- **Major** releases change the benchmark format, execution protocol, or other
  compatibility guarantees.
- **Minor** releases add or remove tasks, or intentionally change evaluation
  behavior in a way that affects score comparability.
- **Patch** releases contain metadata, documentation, infrastructure, or task
  repairs that preserve the intended benchmark definition.

The first planned release is `v0.1.0`. Release tags must point to commits on
`main` and must never be moved after publication.

## One-time Zenodo setup

<!-- TODO(owner): the three `OWNER/longevity-bench` / `harbor-framework`
     references in this section are Terminal-Bench-Science's organisation. This
     repository has no Zenodo deposit and no `CITATION.cff` yet, so treat this
     whole section as a procedure to adapt, not one to run. -->

1. Sign in to [Zenodo](https://zenodo.org/) with a GitHub account that can
   access `OWNER/longevity-bench` (substitute this repository's owner/name — the
   upstream value was `harbor-framework/terminal-bench-science`).
2. If your organisation restricts third-party OAuth apps, an organization
   owner must [approve the Zenodo OAuth
   app](https://docs.github.com/en/organizations/managing-oauth-access-to-your-organizations-data/approving-oauth-apps-for-your-organization).
   Members who are not owners must request approval first.
3. Open the Zenodo
   [GitHub settings](https://zenodo.org/account/settings/github/), synchronize
   repositories, and enable `OWNER/longevity-bench`.
4. Keep [`CITATION.cff`](CITATION.cff) as the single source of release
   metadata. Do not add `.zenodo.json` unless Zenodo-specific grants,
   communities, or relations are needed; if both files exist, Zenodo ignores
   `CITATION.cff`.
   <!-- TODO(owner): CITATION.cff does not exist in this repository yet; the
        upstream file was not ported. Create it and set its concept DOI to your
        own Zenodo record before the archival steps below mean anything. -->

Enabling the repository does not publish anything. Zenodo creates a deposit
only when a GitHub release is published.

## Prepare a release

1. Confirm `main` is current, the working tree is clean, and required CI checks
   pass.
2. Freeze the task set and review all score-affecting changes since the previous
   release.
3. Update [`CITATION.cff`](CITATION.cff) in the release commit:

   ```yaml
   version: v0.1.0
   date-released: YYYY-MM-DD
   ```

4. Validate the citation metadata:

   ```bash
   uvx --from cffconvert cffconvert --validate
   ```

5. Generate the release-note draft:

   ```bash
   uv run --python 3.11 tools/release-notes/generate.py \
     --version v0.1.0 \
     --output release-notes-v0.1.0.md \
     --strict
   ```

   Review the generated task list and add:
   - a short release summary;
   - compatibility or evaluation changes;
   - known issues, if any;
   - citation instructions.
6. Commit the release metadata and notes to `main`, then wait for CI to pass.

## Publish the GitHub release

Create the tag and release from the reviewed commit:

```bash
gh release create v0.1.0 \
  --target main \
  --title "Longevity-Bench v0.1.0" \
  --notes-file release-notes-v0.1.0.md
```

Check that the tag resolves to the intended commit and that the generated
source archives are available. If a published release needs correction, create
a new patch release instead of moving or replacing its tag.

## Verify the Zenodo archive

After GitHub publishes the release:

1. Wait for Zenodo to finish processing the release.
2. Check the title, creators, version, release date, license, description, and
   files on the Zenodo record.
3. Record both identifiers:
   - the **version DOI** identifies that exact benchmark snapshot;
   - the **concept DOI** represents all versions and resolves to the latest one.
4. Confirm the record reports successful archival in Software Heritage.

For the first archived release, add a Zenodo DOI badge to the README using the
concept DOI, add the DOI to the citation guidance, and link the exact version
DOI from the GitHub release notes. Commit those follow-up changes normally; do
not modify the published release tag.
