# frozen_string_literal: true

source "https://rubygems.org"

git_source(:github) {|repo_name| "https://github.com/#{repo_name}" }

gem "github-pages", group: :jekyll_plugins

if Gem::Version.new(RUBY_VERSION) >= Gem::Version.new("3.0.0")
    gem "webrick", ">= 1.6.1"
end
