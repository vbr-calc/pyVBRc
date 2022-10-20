blackv=$(cat .pre-commit-config.yaml | shyaml get-value repos.0.rev | xargs)
isortv=$(cat .pre-commit-config.yaml | shyaml get-value repos.1.rev | xargs)
flk8v=$(cat .pre-commit-config.yaml | shyaml get-value repos.2.rev | xargs)

echo "pre-commit" >> requirements_style.txt
echo "black==$blackv" >> requirements_style.txt
echo "isort==$isortv" >> requirements_style.txt
echo "flake8==$flk8v" >> requirements_style.txt
