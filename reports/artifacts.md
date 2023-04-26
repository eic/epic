---
layout: default
title: Artifacts
permalink: /artifacts
---
{% for artifact in site.static_files %}
  {% if artifact.path contains 'artifacts/' %}
<p><a href="{{ site.baseurl }}/{{ artifact.path }}">{{ artifact.path }}</a></p>
  {% endif %}
{% endfor %}
