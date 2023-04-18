---
layout: default
title: Views (brycecanyon)
permalink: /brycecayon_views
---
{% for image in site.static_files %}
  {% if image.path contains 'epic_brycecanyon_views/' %}
<p><a href="{{ site.baseurl }}/{{ image.path }}">{{ image.path }}</a></p>
    {% if image.extname contains 'png' %}
<img src="{{ site.baseurl }}/{{ image.path }}" width="400px"/>
    {% endif %}
  {% endif %}
{% endfor %}
