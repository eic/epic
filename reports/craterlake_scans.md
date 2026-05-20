---
layout: default
title: Scans (craterlake)
permalink: /craterlake_scans
---
{% for image in site.static_files %}
  {% if image.path contains 'epic_craterlake_scans/' %}
<p><a href="{{ site.baseurl }}/{{ image.path }}">{{ image.path }}</a></p>
    {% if image.extname contains 'png' %}
<img src="{{ site.baseurl }}/{{ image.path }}" width="400px"/>
    {% endif %}
  {% endif %}
{% endfor %}
