if [[ "$(jq '.name' package.json | sed -E 's/(^"|"$)//g')" != "@datagrok/diff-studio-tools" ]]; then
    echo "::error title=ds-tools: failed properties check::Library should be in '@datagrok/diff-studio-tools' scope. Change library name to '@datagrok-libraries/<name>' in /package.json"
fi
