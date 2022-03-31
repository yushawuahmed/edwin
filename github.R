library(usethis)
?use_github
create_github_token()
edit_r_environ()
use_github(protocol = 'https', auth_token = Sys.getenv("GITHUB_PAT"))
